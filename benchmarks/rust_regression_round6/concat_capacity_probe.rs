//! Untimed allocation/capacity replay; no instrumentation in production builds.
use hashrope::{Arena, Node};
use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering::Relaxed};

struct ProbeAllocator;
static ACTIVE: AtomicBool = AtomicBool::new(false);
// alloc, realloc, dealloc, requested bytes, old realloc bytes, new realloc bytes.
static COUNTS: [AtomicU64; 6] = [const { AtomicU64::new(0) }; 6];
static EVENT_COUNT: AtomicUsize = AtomicUsize::new(0);
static EVENTS: [[AtomicU64; 3]; 32] = [const { [const { AtomicU64::new(0) }; 3] }; 32];
unsafe impl GlobalAlloc for ProbeAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if ACTIVE.load(Relaxed) { COUNTS[0].fetch_add(1, Relaxed); COUNTS[3].fetch_add(layout.size() as u64, Relaxed); }
        System.alloc(layout)
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        if ACTIVE.load(Relaxed) { COUNTS[2].fetch_add(1, Relaxed); }
        System.dealloc(ptr, layout)
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let result = System.realloc(ptr, layout, size);
        if ACTIVE.load(Relaxed) {
            COUNTS[1].fetch_add(1, Relaxed);
            COUNTS[3].fetch_add(size as u64, Relaxed);
            COUNTS[4].fetch_add(layout.size() as u64, Relaxed);
            COUNTS[5].fetch_add(size as u64, Relaxed);
            let event = EVENT_COUNT.fetch_add(1, Relaxed);
            if event < EVENTS.len() {
                EVENTS[event][0].store(layout.size() as u64, Relaxed);
                EVENTS[event][1].store(size as u64, Relaxed);
                EVENTS[event][2].store((ptr != result) as u64, Relaxed);
            }
        }
        result
    }
}
#[global_allocator]
static ALLOCATOR: ProbeAllocator = ProbeAllocator;

fn block(label: &str, sample: i32, count: usize, arena: &mut Arena, left: Node, right: Node, expected: u64) {
    let before_len = arena.node_count();
    let before_cap = arena.diagnostic_node_capacity();
    for counter in &COUNTS { counter.store(0, Relaxed); }
    EVENT_COUNT.store(0, Relaxed);
    let mut root = None;
    ACTIVE.store(true, Relaxed);
    for _ in 0..count { root = black_box(arena.concat(black_box(left), black_box(right))); }
    ACTIVE.store(false, Relaxed);
    let counts: [u64; 6] = std::array::from_fn(|i| COUNTS[i].load(Relaxed));
    let event_count = EVENT_COUNT.load(Relaxed);
    assert!(event_count <= EVENTS.len());
    let events: Vec<[u64; 3]> = (0..event_count).map(|i| std::array::from_fn(|j| EVENTS[i][j].load(Relaxed))).collect();
    assert_eq!(arena.node_count(), before_len + count);
    assert_eq!(arena.len(root), 131072);
    assert_eq!(arena.weight(root), 256);
    assert_eq!(arena.hash(root), expected);
    println!("{label}\t{sample}\t{count}\t{before_len}\t{before_cap}\t{}\t{}\t{}\t{}\t{}\t{counts:?}\t{events:?}",
        arena.node_count(), arena.diagnostic_node_capacity(), arena.len(root), arena.weight(root), arena.hash(root));
}

fn main() {
    let data: Vec<u8> = (0..65536).map(|i| (i % 256) as u8).collect();
    let mut joined = data.clone();
    joined.extend_from_slice(&data);
    for (label, warmup) in [("original_growth", 3), ("prewarmed_no_growth", 50003)] {
        let mut arena = Arena::new();
        let left = arena.from_bytes(&data);
        let right = arena.from_bytes(&data);
        assert_eq!(arena.node_count(), 510);
        let expected = arena.hash_bytes(&joined);
        block(label, -1, warmup, &mut arena, left, right, expected);
        for sample in 0..7 { block(label, sample, 2000, &mut arena, left, right, expected); }
    }
}
