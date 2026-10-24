// Diagnostic only: time the same symbolic queries and count query allocations.
use hashrope::{Arena, Node};
use std::{alloc::{GlobalAlloc, Layout, System}, hint::black_box,
    sync::atomic::{AtomicBool, AtomicU64, Ordering}, time::Instant};
struct Alloc;
static COUNT: AtomicBool = AtomicBool::new(false);
static ALLOCS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);
unsafe impl GlobalAlloc for Alloc {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        if COUNT.load(Ordering::Relaxed) { ALLOCS.fetch_add(1,Ordering::Relaxed); BYTES.fetch_add(l.size() as u64,Ordering::Relaxed); }
        System.alloc(l)
    }
    unsafe fn dealloc(&self, p:*mut u8,l:Layout) { System.dealloc(p,l) }
    unsafe fn realloc(&self,p:*mut u8,l:Layout,n:usize)->*mut u8 {
        if COUNT.load(Ordering::Relaxed) { ALLOCS.fetch_add(1,Ordering::Relaxed); BYTES.fetch_add(n as u64,Ordering::Relaxed); }
        System.realloc(p,l,n)
    }
}
#[global_allocator] static GLOBAL: Alloc = Alloc;
fn dag(a:&mut Arena, depth:usize, pattern:&[u8])->Node {
    let mut n=a.from_bytes(pattern);
    for _ in 0..depth { n=a.concat(n,n); }
    n
}
fn ternary(a:&mut Arena, depth:usize)->Node {
    let mut u=a.from_bytes(b"ab"); let mut v=a.from_bytes(b"ba");
    for _ in 0..depth {
        let uv=a.concat(u,v); let vu=a.concat(v,u);
        (u,v)=(a.concat(uv,u),a.concat(vu,v));
    }
    u
}
fn measure(name:&str,a:&Arena,x:Node,y:Node,offset:u64,len:u64,expected:u64) {
    // COUNTER_RESET
    ALLOCS.store(0,Ordering::Relaxed); BYTES.store(0,Ordering::Relaxed);
    COUNT.store(true,Ordering::Relaxed);
    assert_eq!(a.lcp_range(x,offset,y,offset,len),expected);
    COUNT.store(false,Ordering::Relaxed);
    let allocs=ALLOCS.load(Ordering::Relaxed); let bytes=BYTES.load(Ordering::Relaxed);
    let counters=COUNTER_READ;
    let mut times=Vec::new();
    for _ in 0..9 {
        let start=Instant::now();
        for _ in 0..2000 {
            assert_eq!(black_box(a).lcp_range(black_box(x),offset,black_box(y),offset,len),expected);
        }
        times.push(start.elapsed().as_secs_f64()*1e6/2000.0);
    }
    times.sort_by(f64::total_cmp);
    println!("{{\"fixture\":\"{name}\",\"median_us\":{},\"samples_us\":{:?},\"allocations\":{allocs},\"allocated_bytes\":{bytes},\"counters\":{:?}}}",times[4],times,counters);
}
fn main() {
    let mut a=Arena::new();
    let x=dag(&mut a,40,b"ab"); let y=dag(&mut a,40,b"ab");
    measure("doubling_2TiB",&a,x,y,0,1u64<<41,1u64<<41);
    measure("doubling_shifted",&a,x,y,17,(1u64<<41)-33,(1u64<<41)-33);
    let (prefix,_)=a.split(y,(1u64<<41)-1); let end=a.from_bytes(b"?");
    let late=a.concat(prefix,end);
    measure("doubling_late_mismatch",&a,x,late,0,1u64<<41,(1u64<<41)-1);
    let p=a.from_bytes(b"ab"); let q=a.from_bytes(b"abab");
    let x=a.repeat(p,1u64<<40); let y=a.repeat(q,1u64<<39);
    measure("repeat_2TiB",&a,x,y,0,1u64<<41,1u64<<41);
    let x=ternary(&mut a,20); let y=ternary(&mut a,20); let len=2*3u64.pow(20);
    measure("nonperiodic_shared_DAG",&a,x,y,0,len,len);
}
