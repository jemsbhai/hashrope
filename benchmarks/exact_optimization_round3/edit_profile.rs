//! Diagnostic decomposition of the unchanged seeded edit/target workload.
//! Counter builds are separate from timing builds; setup/destruction are untimed.
use hashrope::{Arena, Node, PolynomialHash};
use ropey::Rope;
use std::{hint::black_box, time::Instant};

#[cfg(feature = "counters")]
mod allocations {
    use std::alloc::{GlobalAlloc, Layout, System};
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering::Relaxed};
    pub struct Allocator;
    static ACTIVE: AtomicBool = AtomicBool::new(false);
    static VALUES: [AtomicU64; 5] = [const { AtomicU64::new(0) }; 5];
    unsafe impl GlobalAlloc for Allocator {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            let result = System.alloc(layout);
            if ACTIVE.load(Relaxed) && !result.is_null() {
                VALUES[0].fetch_add(1, Relaxed);
                VALUES[1].fetch_add(layout.size() as u64, Relaxed);
            }
            result
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            let result = System.alloc_zeroed(layout);
            if ACTIVE.load(Relaxed) && !result.is_null() {
                VALUES[0].fetch_add(1, Relaxed);
                VALUES[1].fetch_add(layout.size() as u64, Relaxed);
            }
            result
        }
        unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, size: usize) -> *mut u8 {
            let result = System.realloc(pointer, layout, size);
            if ACTIVE.load(Relaxed) && !result.is_null() {
                VALUES[2].fetch_add(1, Relaxed);
                VALUES[3].fetch_add(size as u64, Relaxed);
            }
            result
        }
        unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
            if ACTIVE.load(Relaxed) { VALUES[4].fetch_add(1, Relaxed); }
            System.dealloc(pointer, layout)
        }
    }
    pub fn start() {
        for value in &VALUES { value.store(0, Relaxed); }
        ACTIVE.store(true, Relaxed);
    }
    pub fn end() -> [u64; 5] {
        ACTIVE.store(false, Relaxed);
        std::array::from_fn(|i| VALUES[i].load(Relaxed))
    }
}
#[cfg(feature = "counters")]
#[global_allocator]
static ALLOCATOR: allocations::Allocator = allocations::Allocator;

struct Measurement { nanoseconds: u128, allocations: [u64; 5], hashes: [u64; 9], answer: u64 }
fn measure(mut operation: impl FnMut() -> u64) -> Measurement {
    #[cfg(feature = "counters")]
    { hashrope::diagnostic::start(); allocations::start(); }
    let started = Instant::now();
    let answer = black_box(operation());
    let nanoseconds = started.elapsed().as_nanos();
    #[cfg(feature = "counters")]
    let (allocations, hashes) = (allocations::end(), hashrope::diagnostic::end());
    #[cfg(not(feature = "counters"))]
    let (allocations, hashes) = ([0; 5], [0; 9]);
    Measurement { nanoseconds, allocations, hashes, answer }
}
fn report(size: usize, shared: bool, arm: &str, phase: &str, sample: usize,
          operations: usize, nodes: usize, value: Measurement) {
    print!("{size}\t{shared}\t{arm}\t{phase}\t{sample}\t{operations}\t{}\t{nodes}", value.nanoseconds);
    for count in value.allocations { print!("\t{count}"); }
    for count in value.hashes { print!("\t{count}"); }
    println!();
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize { self.next() as usize % n }
    fn ascii(&mut self) -> u8 { 32 + self.below(95) as u8 }
}
struct Edit { delete: usize, insert: usize, bytes: Vec<u8> }
struct Workload { initial: Vec<u8>, edits: Vec<Edit>, states: Vec<Vec<u8>>, targets: Vec<Vec<u8>> }
fn workload(n: usize, seed: u64) -> Workload {
    let mut rng = Rng(seed);
    let initial: Vec<u8> = (0..n).map(|_| rng.ascii()).collect();
    let mut oracle = initial.clone();
    let (mut edits, mut states, mut targets) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..50 {
        let len = 1 + rng.below(32.min(n - 1));
        let edit = Edit { delete: rng.below(n - len + 1), insert: rng.below(n - len + 1),
                          bytes: (0..len).map(|_| rng.ascii()).collect() };
        oracle.drain(edit.delete..edit.delete + len);
        oracle.splice(edit.insert..edit.insert, edit.bytes.iter().copied());
        states.push(oracle.clone());
        let mut target = oracle.clone();
        if i % 2 == 1 { target[n - 1] = if target[n - 1] == b'!' { b'?' } else { b'!' }; }
        edits.push(edit); targets.push(target);
    }
    Workload { initial, edits, states, targets }
}
fn edit_hashrope(arena: &mut Arena, root: Node, edit: &Edit) -> Node {
    let (left, rest) = arena.split(root, edit.delete as u64);
    let (_, right) = arena.split(rest, edit.bytes.len() as u64);
    let root = arena.concat(left, right);
    let (left, right) = arena.split(root, edit.insert as u64);
    let added = arena.from_bytes(&edit.bytes);
    let head = arena.concat(left, added);
    arena.concat(head, right)
}
fn edit_ropey(rope: &mut Rope, edit: &Edit) {
    rope.remove(edit.delete..edit.delete + edit.bytes.len());
    rope.insert(edit.insert, std::str::from_utf8(&edit.bytes).unwrap());
}
fn equal_mask() -> u64 { (0..50).step_by(2).map(|i| 1u64 << i).sum() }

fn hashrope_case(w: &Workload, shared: bool, lazy: bool, phase: &str, sample: usize) {
    let create = || if lazy { Arena::new_lazy() } else { Arena::new() };
    let (mut arena, mut independent) = (create(), create());
    let mut root = arena.from_bytes(&w.initial);
    let mut targets = Vec::new();
    if shared {
        let mut target = root;
        for (i, edit) in w.edits.iter().enumerate() {
            target = edit_hashrope(&mut arena, target, edit);
            let copied = if i % 2 == 1 {
                let (prefix, _) = arena.split(target, w.initial.len() as u64 - 1);
                let ending = arena.from_bytes(&w.targets[i][w.initial.len() - 1..]);
                arena.concat(prefix, ending)
            } else { target };
            targets.push(copied);
        }
    } else {
        for target in &w.targets { targets.push(independent.from_bytes(target)); }
    }
    // Identical to the existing harness's untimed Hashrope edit warmup.
    let mut warm = root;
    for edit in &w.edits { warm = edit_hashrope(&mut arena, warm, edit); }
    black_box(warm);
    let mut states = Vec::new();
    if phase == "compare_only" || phase == "fingerprint_only" {
        for edit in &w.edits {
            root = edit_hashrope(&mut arena, root, edit);
            states.push(root);
        }
    }
    let before = arena.node_count();
    let repeats = if phase == "compare_only" || phase == "fingerprint_only" { 3 } else { 1 };
    let measured = measure(|| {
        let mut all_masks = 0;
        for _ in 0..repeats {
            let mut mask = 0;
            for (i, edit) in w.edits.iter().enumerate() {
                if phase == "edit_only" || phase == "combined" { root = edit_hashrope(&mut arena, root, edit); }
                if phase != "edit_only" {
                    let state = if phase == "combined" { root } else { states[i] };
                    let equal = if phase == "fingerprint_only" {
                        arena.hash(state) == if shared { arena.hash(targets[i]) } else { independent.hash(targets[i]) }
                    } else if shared { arena.equal(state, targets[i]) }
                    else { arena.equal_across(state, &independent, targets[i]) };
                    mask |= (black_box(equal) as u64) << i;
                }
            }
            all_masks ^= black_box(mask);
        }
        all_masks
    });
    if phase != "edit_only" { assert_eq!(measured.answer, equal_mask()); }
    assert_eq!(arena.to_bytes(root), *w.states.last().unwrap());
    for (i, state) in states.iter().enumerate() { assert_eq!(arena.to_bytes(*state), w.states[i]); }
    report(w.initial.len(), shared, if lazy { "hashrope_lazy" } else { "hashrope_eager" }, phase,
           sample, 50 * repeats, arena.node_count() - before, measured);
}

fn ropey_case(w: &Workload, shared: bool, phase: &str, sample: usize) {
    let mut rope = Rope::from_str(std::str::from_utf8(&w.initial).unwrap());
    let mut targets = Vec::new();
    if shared {
        let mut target = rope.clone();
        for (i, edit) in w.edits.iter().enumerate() {
            edit_ropey(&mut target, edit);
            let mut copy = target.clone();
            if i % 2 == 1 {
                let n = w.initial.len();
                copy.remove(n - 1..n);
                copy.insert(n - 1, std::str::from_utf8(&w.targets[i][n - 1..]).unwrap());
            }
            targets.push(copy);
        }
    } else {
        for target in &w.targets { targets.push(Rope::from_str(std::str::from_utf8(target).unwrap())); }
    }
    let mut states = Vec::new();
    if phase == "compare_only" {
        for edit in &w.edits { edit_ropey(&mut rope, edit); states.push(rope.clone()); }
    }
    let repeats = if phase == "compare_only" { 3 } else { 1 };
    let measured = measure(|| {
        let mut all_masks = 0;
        for _ in 0..repeats {
            let mut mask = 0;
            for (i, edit) in w.edits.iter().enumerate() {
                if phase != "compare_only" { edit_ropey(&mut rope, edit); }
                if phase != "edit_only" {
                    let state = if phase == "combined" { &rope } else { &states[i] };
                    mask |= (black_box(state == &targets[i]) as u64) << i;
                }
            }
            all_masks ^= black_box(mask);
        }
        all_masks
    });
    if phase != "edit_only" { assert_eq!(measured.answer, equal_mask()); }
    assert_eq!(rope.to_string().as_bytes(), w.states.last().unwrap());
    for (i, state) in states.iter().enumerate() { assert_eq!(state.to_string().as_bytes(), w.states[i]); }
    report(w.initial.len(), shared, "ropey", phase, sample, 50 * repeats, 0, measured);
}

fn primitives(samples: usize) {
    let mut hasher = PolynomialHash::default_hash();
    for length in [1usize, 16, 256, 512] {
        let bytes: Vec<u8> = (0..length).map(|i| (i * 73) as u8).collect();
        let expected = hasher.hash(&bytes);
        let count = if cfg!(feature = "counters") { 1 } else { 100_000 };
        for sample in 0..samples {
            let measured = measure(|| {
                let mut value = 0;
                for _ in 0..count { value = hasher.hash(black_box(&bytes)); black_box(value); }
                value
            });
            assert_eq!(measured.answer, expected);
            report(length, false, "hashrope_eager", "primitive_hash_bytes", sample, count, 0, measured);
        }
        let expected = hasher.hash_concat(123, length as u64, 456);
        for sample in 0..samples {
            let measured = measure(|| {
                let mut value = 0;
                for _ in 0..count { value = hasher.hash_concat(black_box(123), black_box(length as u64), black_box(456)); black_box(value); }
                value
            });
            assert_eq!(measured.answer, expected);
            report(length, false, "hashrope_eager", "primitive_hash_concat_cached", sample, count, 0, measured);
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seed: u64 = args[1].parse().unwrap();
    let samples: usize = args.get(2).map(|x| x.parse().unwrap()).unwrap_or(3);
    primitives(samples);
    let sizes: Vec<usize> = args.get(3).map(|value| value.split(',').map(|n| n.parse().unwrap()).collect())
        .unwrap_or_else(|| vec![10_000, 100_000, 1_000_000]);
    for size in sizes {
        let work = workload(size, seed);
        for shared in [false, true] {
            let mut arms = ["hashrope_eager", "hashrope_lazy", "ropey"];
            arms.rotate_left(seed as usize % 3);
            for sample in 0..samples {
                let mut phases = ["edit_only", "compare_only", "combined"];
                phases.rotate_left((seed as usize + sample) % 3);
                for phase in phases {
                    for arm in arms {
                        if arm == "ropey" { ropey_case(&work, shared, phase, sample); }
                        else { hashrope_case(&work, shared, arm == "hashrope_lazy", phase, sample); }
                    }
                }
                hashrope_case(&work, shared, false, "fingerprint_only", sample);
            }
        }
    }
}
