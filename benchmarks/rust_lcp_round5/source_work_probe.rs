//! Source counters and requested allocator traffic only; no timed measurements.
use hashrope::{Arena, Node};
use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};

struct CountAlloc;
static ACTIVE: AtomicBool = AtomicBool::new(false);
static ALLOC: AtomicU64 = AtomicU64::new(0);
static REALLOC: AtomicU64 = AtomicU64::new(0);
static DEALLOC: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);
unsafe impl GlobalAlloc for CountAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if ACTIVE.load(Relaxed) { ALLOC.fetch_add(1, Relaxed); BYTES.fetch_add(layout.size() as u64, Relaxed); }
        System.alloc(layout)
    }
    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        if ACTIVE.load(Relaxed) { DEALLOC.fetch_add(1, Relaxed); }
        System.dealloc(p, layout)
    }
    unsafe fn realloc(&self, p: *mut u8, layout: Layout, new: usize) -> *mut u8 {
        if ACTIVE.load(Relaxed) { REALLOC.fetch_add(1, Relaxed); BYTES.fetch_add(new as u64, Relaxed); }
        System.realloc(p, layout, new)
    }
}
#[global_allocator]
static GLOBAL: CountAlloc = CountAlloc;

fn probe(name: &str, arena: &Arena, a: Node, ao: u64, b: Node, bo: u64, length: u64, expected: u64) {
    hashrope::rope::source_probe_reset();
    for value in [&ALLOC, &REALLOC, &DEALLOC, &BYTES] { value.store(0, Relaxed); }
    ACTIVE.store(true, Relaxed);
    let answer = black_box(arena.lcp_range(black_box(a), ao, black_box(b), bo, length));
    ACTIVE.store(false, Relaxed);
    let counts = hashrope::rope::source_probe_counts();
    let alloc = ALLOC.load(Relaxed);
    let realloc = REALLOC.load(Relaxed);
    let dealloc = DEALLOC.load(Relaxed);
    let bytes = BYTES.load(Relaxed);
    assert_eq!(answer, expected, "{name}");
    println!("{name}\t{answer}\t{expected}\t{alloc}\t{realloc}\t{dealloc}\t{bytes}\t{counts:?}");
}

fn main() {
    for n in [65536usize, 1048576] {
        let bytes: Vec<_> = (0..n).map(|x| (x % 256) as u8).collect();
        let mut arena = Arena::new();
        let a = arena.from_bytes(&bytes);
        let b = arena.from_bytes(&bytes);
        probe(&format!("aligned/{n}"), &arena, a, 0, b, 0, n as u64, n as u64);
        let mut changed = bytes;
        changed[n - 1] ^= 1;
        let b = arena.from_bytes(&changed);
        probe(&format!("aligned_late/{n}"), &arena, a, 0, b, 0, n as u64, n as u64 - 1);
    }
    let mut arena = Arena::new();
    let mut a = arena.from_bytes(b"ab");
    let mut b = arena.from_bytes(b"ab");
    for _ in 0..40 { a = arena.concat(a, a); b = arena.concat(b, b); }
    probe("self_concat_2TiB", &arena, a, 0, b, 0, 1u64 << 41, 1u64 << 41);
    probe("self_concat_shifted", &arena, a, 5, b, 7, (1u64 << 41) - 7, (1u64 << 41) - 7);

    let mut arena = Arena::new();
    let a_byte = arena.from_bytes(b"a");
    let b_byte = arena.from_bytes(b"b");
    let n = 1u64 << 32;
    let d = 17u64;
    let aa = arena.repeat(a_byte, n);
    let bb = arena.repeat(b_byte, n);
    let unit_a = arena.concat(aa, bb);
    let independent_a = arena.from_bytes(b"a");
    let independent_b = arena.from_bytes(b"b");
    let short_a = arena.repeat(independent_a, n - d);
    let tiny_a = arena.repeat(independent_a, d);
    let independent_bb = arena.repeat(independent_b, n);
    let bb_tiny = arena.concat(independent_bb, tiny_a);
    let unit_b = arena.concat(short_a, bb_tiny);
    let a = arena.repeat(unit_a, 8);
    let b = arena.repeat(unit_b, 8);
    // Independent units rotated by17 bytes, both phases nonzero. One complete
    // period is partitioned at the two different wrap points, giving3 pieces.
    probe("compressed_unit_three_pieces", &arena, a, d + 1, b, 1, 4 * n, 4 * n);
    let (_, tail) = arena.split(unit_b, 1);
    let c_byte = arena.from_bytes(b"c");
    let changed_unit = arena.concat(c_byte, tail);
    let changed = arena.repeat(changed_unit, 8);
    probe("compressed_unit_last_witness_byte", &arena, a, d + 1, changed, 1, 4 * n, 2 * n - 1);
}
