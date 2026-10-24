//! Separate controlled diagnostic; no instrumentation in either timed core.
use std::{hint::black_box, time::{Duration, Instant}};

#[cfg(windows)]
mod platform {
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> *mut std::ffi::c_void;
        fn GetProcessAffinityMask(p: *mut std::ffi::c_void, process: *mut usize, system: *mut usize) -> i32;
        fn SetProcessAffinityMask(p: *mut std::ffi::c_void, mask: usize) -> i32;
        fn GetCurrentProcessorNumber() -> u32;
    }
    pub fn processor() -> u32 { unsafe { GetCurrentProcessorNumber() } }
    pub fn affinity(requested: usize) -> (usize, usize) {
        unsafe {
            let p = GetCurrentProcess();
            let (mut initial, mut system) = (0, 0);
            assert_ne!(GetProcessAffinityMask(p, &mut initial, &mut system), 0);
            if requested != 0 {
                assert_eq!(requested & initial, requested);
                assert_ne!(SetProcessAffinityMask(p, requested), 0);
            }
            let mut actual = 0;
            assert_ne!(GetProcessAffinityMask(p, &mut actual, &mut system), 0);
            (initial, actual)
        }
    }
}
#[cfg(not(windows))]
mod platform {
    pub fn processor() -> u32 { u32::MAX }
    pub fn affinity(requested: usize) -> (usize, usize) { assert_eq!(requested, 0); (0, 0) }
}

#[derive(Clone, Copy, Debug)]
enum Case { ConcatGrowing, LcpIndependent, Combined, EditOnly, CompareOnly }
impl Case {
    fn name(self) -> &'static str {
        match self {
            Self::ConcatGrowing => "concat_growing", Self::LcpIndependent => "lcp_independent_1MiB",
            Self::Combined => "churn_combined_100K", Self::EditOnly => "churn_edit_only_100K",
            Self::CompareOnly => "churn_compare_only_100K",
        }
    }
    fn index(self) -> usize { self as usize }
    fn concat(self) -> bool { matches!(self, Self::ConcatGrowing) }
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
struct Raw { initial: Vec<u8>, operations: Vec<Edit>, targets: Vec<Vec<u8>>, final_bytes: Vec<u8>, concat: Vec<u8>, lcp: Vec<u8> }
impl Raw {
    fn new(seed: u64) -> Self {
        let n = 100_000;
        let mut rng = Rng(seed);
        let initial: Vec<u8> = (0..n).map(|_| rng.ascii()).collect();
        let mut oracle = initial.clone();
        let mut operations = Vec::new();
        let mut targets = Vec::new();
        for i in 0..50 {
            let len = 1 + rng.below(32.min(n - 1));
            let edit = Edit { delete: rng.below(n - len + 1), insert: rng.below(n - len + 1),
                bytes: (0..len).map(|_| rng.ascii()).collect() };
            oracle.drain(edit.delete..edit.delete + len);
            oracle.splice(edit.insert..edit.insert, edit.bytes.iter().copied());
            let mut target = oracle.clone();
            if i % 2 == 1 { target[n - 1] = if target[n - 1] == b'!' { b'?' } else { b'!' }; }
            operations.push(edit);
            targets.push(target);
        }
        Self { initial, operations, targets, final_bytes: oracle,
               concat: (0..65536).map(|i| (i % 256) as u8).collect(),
               lcp: (0..1_048_576).map(|i| (i % 256) as u8).collect() }
    }
}
fn expected_mask() -> u64 { (0..50).step_by(2).map(|i| 1u64 << i).sum() }
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    if sorted.len() % 2 == 1 { sorted[sorted.len() / 2] }
    else { (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0 }
}
struct Session {
    durations: Vec<u128>, operations_per_sample: usize, answer: u64,
    nodes_start: usize, nodes_end: usize, root_address: usize, leaf_address: usize,
}
impl Session {
    fn statistic(&self) -> f64 {
        median(&self.durations.iter().map(|n| *n as f64 / self.operations_per_sample as f64).collect::<Vec<_>>())
    }
}
struct Block { sessions: Vec<Session>, wall_ns: u128, cpu_before: u32, cpu_after: u32 }
impl Block {
    fn timed_ns(&self) -> u128 { self.sessions.iter().flat_map(|s| &s.durations).sum() }
    fn statistic(&self) -> f64 { median(&self.sessions.iter().map(Session::statistic).collect::<Vec<_>>()) }
    fn operations(&self) -> usize { self.sessions.iter().map(|s| s.operations_per_sample * s.durations.len()).sum() }
    fn emit(&self, case: Case, order: usize, position: usize, slot: &str) {
        for (index, session) in self.sessions.iter().enumerate() {
            for (sample, ns) in session.durations.iter().enumerate() {
                println!("sample\t{}\t{order}\t{position}\t{slot}\t{index}\t{sample}\t{}\t{ns}\t{}\t{}\t{}\t{}\t{}",
                    case.name(), session.operations_per_sample, session.answer,
                    session.nodes_start, session.nodes_end, session.root_address, session.leaf_address);
            }
        }
        println!("block\t{}\t{order}\t{position}\t{slot}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            case.name(), self.sessions.len(), self.operations(), self.statistic(), self.timed_ns(),
            self.wall_ns, self.cpu_before, self.cpu_after);
    }
}

macro_rules! implementation {
    ($module:ident, $core:ident) => {
        mod $module {
            use super::*;
            use $core::{Arena, Node, NodeInner};
            fn edit(a: &mut Arena, root: Node, operation: &Edit) -> Node {
                let (left, rest) = a.split(root, operation.delete as u64);
                let (_, right) = a.split(rest, operation.bytes.len() as u64);
                let root = a.concat(left, right);
                let (left, right) = a.split(root, operation.insert as u64);
                let added = a.from_bytes(&operation.bytes);
                let head = a.concat(left, added);
                a.concat(head, right)
            }
            fn addresses(arena: &Arena, root: Node) -> (usize, usize) {
                let mut id = root.unwrap();
                let root_address = arena.node(id) as *const NodeInner as usize;
                loop {
                    match arena.node(id) {
                        NodeInner::Leaf { data, .. } => return (root_address, data.as_ptr() as usize),
                        NodeInner::Internal { left, .. } => id = *left,
                        NodeInner::Repeat { child, .. } => id = *child,
                    }
                }
            }
            pub fn collision_oracle() {
                let mut a = Arena::new();
                let x = a.from_bytes(&[0, 131]);
                let y = a.from_bytes(&[1, 0]);
                assert_eq!(a.hash(x), a.hash(y));
                assert!(!a.equal(x, y));
                assert_eq!(a.lcp(x, y), 0);
            }
            fn hash_lcp(arena: &mut Arena, a: Node, b: Node) -> u64 {
                let (mut lo, mut hi) = (0, arena.len(a).min(arena.len(b)));
                while lo < hi {
                    let mid = lo + (hi - lo + 1) / 2;
                    if arena.substr_hash(a, 0, mid) == arena.substr_hash(b, 0, mid) {
                        lo = mid;
                    } else { hi = mid - 1; }
                }
                lo
            }
            fn lcp_session(raw: &Raw) -> Session {
                // Match the original microbenchmark's complete local arena
                // construction, rather than allocating only the two equal roots.
                let size = raw.lcp.len();
                let mut arena = Arena::new();
                let first = arena.from_bytes(&raw.lcp);
                let independent = arena.from_bytes(&raw.lcp);
                let mut changed = raw.lcp.clone();
                changed[0] ^= 1;
                let early = arena.from_bytes(&changed);
                changed[0] ^= 1;
                changed[size - 1] ^= 1;
                let late = arena.from_bytes(&changed);
                let (prefix, _) = arena.split(first, size as u64 - 1);
                let ending = arena.from_bytes(&changed[size - 1..]);
                let shared_late = arena.concat(prefix, ending);
                assert_ne!(first, independent);
                assert_ne!(addresses(&arena, first).1, addresses(&arena, independent).1);
                assert_eq!(arena.to_bytes(first), raw.lcp);
                assert_eq!(arena.to_bytes(independent), raw.lcp);
                let early_bytes = arena.to_bytes(early);
                assert_eq!(early_bytes[0], raw.lcp[0] ^ 1);
                assert_eq!(early_bytes[1..], raw.lcp[1..]);
                drop(early_bytes);
                assert_eq!(arena.to_bytes(late), changed);
                assert_eq!(arena.to_bytes(shared_late), changed);
                assert_eq!(arena.lcp(first, early), 0);
                assert_eq!(arena.lcp(first, late), size as u64 - 1);
                assert_eq!(arena.lcp(first, shared_late), size as u64 - 1);
                // Preserve the target case's original untimed lcp/equal/hash
                // checks. The preceding same-root timed cells are not replayed.
                assert_eq!(arena.lcp(first, independent), size as u64);
                assert!(arena.equal(first, independent));
                assert_eq!(hash_lcp(&mut arena, first, independent), size as u64);
                for _ in 0..3 {
                    black_box(arena.lcp(black_box(first), black_box(independent)));
                }
                let nodes_start = arena.node_count();
                let (root_address, leaf_address) = addresses(&arena, first);
                let mut durations = Vec::with_capacity(7);
                let mut answers = [0u64; 7];
                for answer in &mut answers {
                    let start = Instant::now();
                    for _ in 0..100 {
                        *answer = black_box(black_box(&arena).lcp(black_box(first), black_box(independent)));
                    }
                    durations.push(start.elapsed().as_nanos());
                }
                assert!(answers.iter().all(|answer| *answer == size as u64));
                assert_eq!(arena.node_count(), nodes_start);
                assert_eq!(arena.to_bytes(first), raw.lcp);
                assert_eq!(arena.to_bytes(independent), raw.lcp);
                Session { durations, operations_per_sample: 100, answer: size as u64,
                          nodes_start, nodes_end: arena.node_count(), root_address, leaf_address }
            }
            fn concat_session(raw: &Raw) -> Session {
                let mut arena = Arena::new();
                let left = arena.from_bytes(&raw.concat);
                let right = arena.from_bytes(&raw.concat);
                assert_eq!(arena.node_count(), 510);
                let warm = 3;
                for _ in 0..warm { black_box(arena.concat(black_box(left), black_box(right))); }
                let nodes_start = arena.node_count();
                assert_eq!(nodes_start, 510 + warm);
                let (root_address, leaf_address) = addresses(&arena, left);
                let mut durations = Vec::with_capacity(7);
                let mut endings = [None; 7];
                for ending in &mut endings {
                    let start = Instant::now();
                    for _ in 0..2000 {
                        *ending = black_box(arena.concat(black_box(left), black_box(right)));
                    }
                    durations.push(start.elapsed().as_nanos());
                }
                let nodes_end = arena.node_count();
                assert_eq!(nodes_end, nodes_start + 14_000);
                for root in endings {
                    assert_eq!(arena.len(root), 131_072);
                    assert_eq!(arena.weight(root), 256);
                    assert_eq!(arena.hash(root), 2_138_634_010_386_205_116);
                    let bytes = arena.to_bytes(root);
                    assert_eq!(&bytes[..65536], raw.concat.as_slice());
                    assert_eq!(&bytes[65536..], raw.concat.as_slice());
                }
                Session { durations, operations_per_sample: 2000, answer: 131_072,
                          nodes_start, nodes_end, root_address, leaf_address }
            }
            fn churn_session(raw: &Raw, case: Case) -> Session {
                let mut arena = Arena::new();
                let mut root = arena.from_bytes(&raw.initial);
                assert_eq!(arena.to_bytes(root), raw.initial);
                let mut independent = Arena::new();
                let targets: Vec<Node> = if matches!(case, Case::EditOnly) { Vec::new() }
                    else { raw.targets.iter().map(|bytes| independent.from_bytes(bytes)).collect() };
                for (&target, bytes) in targets.iter().zip(&raw.targets) {
                    assert_eq!(independent.to_bytes(target), *bytes);
                }
                let mut warm = root;
                for operation in &raw.operations { warm = edit(&mut arena, warm, operation); }
                black_box(warm);
                let mut states = [None; 50];
                if matches!(case, Case::CompareOnly) {
                    for (i, operation) in raw.operations.iter().enumerate() {
                        root = edit(&mut arena, root, operation);
                        states[i] = root;
                        let bytes = arena.to_bytes(root);
                        assert_eq!(bytes == raw.targets[i], i % 2 == 0);
                        if i % 2 == 1 {
                            assert_eq!(bytes[..99_999], raw.targets[i][..99_999]);
                            assert_ne!(bytes[99_999], raw.targets[i][99_999]);
                        }
                    }
                    assert_eq!(arena.to_bytes(root), raw.final_bytes);
                }
                let nodes_start = arena.node_count();
                let (root_address, leaf_address) = addresses(&arena, root);
                let mut answer = 0u64;
                let start = Instant::now();
                match case {
                    Case::Combined => for (i, operation) in raw.operations.iter().enumerate() {
                        root = edit(&mut arena, black_box(root), operation);
                        let equal = black_box(black_box(&arena).equal_across(black_box(root), black_box(&independent), black_box(targets[i])));
                        answer |= (equal as u64) << i;
                    },
                    Case::EditOnly => for operation in &raw.operations {
                        root = black_box(edit(&mut arena, black_box(root), operation));
                    },
                    Case::CompareOnly => for i in 0..50 {
                        let equal = black_box(black_box(&arena).equal_across(black_box(states[i]), black_box(&independent), black_box(targets[i])));
                        answer |= (equal as u64) << i;
                    },
                    _ => unreachable!(),
                }
                let elapsed = start.elapsed().as_nanos();
                assert_eq!(arena.to_bytes(root), raw.final_bytes);
                if matches!(case, Case::EditOnly) { answer = arena.len(root); }
                else { assert_eq!(answer, expected_mask()); assert_eq!(answer.count_ones(), 25); }
                Session { durations: vec![elapsed], operations_per_sample: 50, answer,
                          nodes_start, nodes_end: arena.node_count(), root_address, leaf_address }
            }
            #[inline(never)]
            pub fn measure(raw: &Raw, case: Case, count: usize) -> Block {
                let cpu_before = platform::processor();
                let start = Instant::now();
                let sessions = (0..count).map(|_| if case.concat() { concat_session(raw) }
                    else if matches!(case, Case::LcpIndependent) { lcp_session(raw) }
                    else { churn_session(raw, case) }).collect();
                let wall_ns = start.elapsed().as_nanos();
                let cpu_after = platform::processor();
                Block { sessions, wall_ns, cpu_before, cpu_after }
            }
        }
    }
}
implementation!(a, slot_a);
implementation!(b, slot_b);

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seed: usize = args[1].parse().unwrap();
    let requested: usize = args[2].parse().unwrap();
    assert!((1..=6).contains(&seed));
    let (initial_affinity, actual_affinity) = platform::affinity(requested);
    let churn_seed = ((seed - 1) % 3 + 1) as u64;
    println!("meta\t{seed}\t{churn_seed}\t{requested}\t{initial_affinity}\t{actual_affinity}");
    a::collision_oracle(); b::collision_oracle();
    let raw = Raw::new(churn_seed);
    let mut cases = [Case::ConcatGrowing, Case::LcpIndependent, Case::Combined, Case::EditOnly, Case::CompareOnly];
    let mut rng = Rng(seed as u64 ^ 0xd103_ab74_8923_71f0);
    for i in (1..cases.len()).rev() { cases.swap(i, rng.below(i + 1)); }
    for (order, case) in cases.into_iter().enumerate() {
        // Independent of shuffled case position and churn seed identity.
        let first_a = (seed + case.index()) % 2 == 0;
        let mut count = 1usize;
        let duration = loop {
            let (left, right) = if first_a {
                (a::measure(&raw, case, count), b::measure(&raw, case, count))
            } else {
                let right = b::measure(&raw, case, count);
                let left = a::measure(&raw, case, count);
                (left, right)
            };
            let elapsed = left.timed_ns().max(right.timed_ns());
            if elapsed >= 1_000_000 { break elapsed; }
            count = count.checked_mul(2).unwrap();
            assert!(count <= 4096, "calibration runaway");
        };
        count = ((count as u128 * 25_000_000) / duration).max(1) as usize;
        assert!(count <= 4096, "calibration budget");
        let warm = Instant::now();
        let mut warm_a_first = first_a;
        while warm.elapsed() < Duration::from_millis(20) {
            if warm_a_first {
                black_box(a::measure(&raw, case, count)); black_box(b::measure(&raw, case, count));
            } else {
                black_box(b::measure(&raw, case, count)); black_box(a::measure(&raw, case, count));
            }
            warm_a_first = !warm_a_first;
        }
        println!("calibration\t{}\t{order}\t{count}\t{}", case.name(), if first_a { "ABBA" } else { "BAAB" });
        let slots = if first_a { [true, false, false, true] } else { [false, true, true, false] };
        for (position, is_a) in slots.into_iter().enumerate() {
            let block = if is_a { a::measure(&raw, case, count) } else { b::measure(&raw, case, count) };
            block.emit(case, order, position, if is_a { "a" } else { "b" });
        }
    }
}
