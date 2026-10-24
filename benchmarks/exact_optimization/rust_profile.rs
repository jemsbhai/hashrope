use hashrope::Arena;
use std::{alloc::{GlobalAlloc, Layout, System}, sync::atomic::{AtomicBool, AtomicU64, Ordering}, hint::black_box, time::Instant};
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
fn main() {
    for n in [10_000,100_000,1_000_000] {
        let bytes: Vec<u8> = (0..n).map(|i| b' '+(i%95) as u8).collect();
        let mut a=Arena::new(); let mut b=Arena::new();
        let x=a.from_bytes(&bytes); let y=b.from_bytes(&bytes);
        // One equality call counts allocation and structural work after setup.
        hashrope::rope::profile_reset();
        ALLOCS.store(0,Ordering::Relaxed); BYTES.store(0,Ordering::Relaxed); COUNT.store(true,Ordering::Relaxed);
        assert!(a.equal_across(x,&b,y));
        COUNT.store(false,Ordering::Relaxed);
        let counts=hashrope::rope::profile_read();
        let allocs=ALLOCS.load(Ordering::Relaxed); let allocated=BYTES.load(Ordering::Relaxed);
        let mut samples=Vec::new();
        for _ in 0..7 {
            let t=Instant::now();
            for _ in 0..100 { assert!(black_box(a.equal_across(x,&b,y))); }
            samples.push(t.elapsed().as_secs_f64()*1e6/100.0);
        }
        samples.sort_by(f64::total_cmp);
        let flat=bytes.clone();
        let t=Instant::now();
        for _ in 0..1000 { assert!(black_box(bytes.as_slice()) == black_box(flat.as_slice())); }
        let flat_us=t.elapsed().as_secs_f64()*1e6/1000.0;
        println!("{{\"bytes\":{n},\"query_allocations\":{allocs},\"requested_allocation_bytes\":{allocated},\"median_us\":{},\"samples_us\":{:?},\"flat_slice_us\":{flat_us},\"counters\":{:?},\"instrumented\":{}}}",samples[3],samples,counts,PROFILE_COUNTERS);
    }
}
