//! Diagnostic only. The original comparison_safety benchmark remains unchanged.
use std::{hint::black_box, time::{Duration, Instant}};

#[cfg(windows)]
mod platform {
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> *mut std::ffi::c_void;
        fn GetProcessAffinityMask(process: *mut std::ffi::c_void, process_mask: *mut usize, system_mask: *mut usize) -> i32;
        fn SetProcessAffinityMask(process: *mut std::ffi::c_void, mask: usize) -> i32;
        fn GetCurrentProcessorNumber() -> u32;
    }
    pub fn processor() -> u32 { unsafe { GetCurrentProcessorNumber() } }
    pub fn affinity(requested: usize) -> (usize, usize) {
        unsafe {
            let process = GetCurrentProcess();
            let (mut initial, mut system) = (0, 0);
            assert_ne!(GetProcessAffinityMask(process, &mut initial, &mut system), 0);
            if requested != 0 { assert_eq!(requested & initial, requested); assert_ne!(SetProcessAffinityMask(process, requested), 0); }
            let mut actual = 0;
            assert_ne!(GetProcessAffinityMask(process, &mut actual, &mut system), 0);
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
enum Operation { Lcp, Equal }
impl Operation { fn name(self) -> &'static str { match self { Self::Lcp => "lcp", Self::Equal => "equal" } } }
struct Row { elapsed_ns: u128, iterations: usize, answer: u64, processor_before: u32, processor_after: u32 }

macro_rules! implementation {
    ($module:ident, $core:ident) => {
        mod $module {
            use super::*;
            use $core::{Arena, Node, NodeInner};
            struct Fixture { arena: usize, first: Node, second: Node, expected: u64, length: u64, name: String }
            pub struct Suite { arenas: Vec<Arena>, fixtures: Vec<Fixture> }
            impl Suite {
                pub fn new() -> Self {
                    let mut result = Self { arenas: Vec::new(), fixtures: Vec::new() };
                    for size in [65536usize, 1048576] {
                        let raw: Vec<u8> = (0..size).map(|i| (i % 256) as u8).collect();
                        let mut arena = Arena::new();
                        let first = arena.from_bytes(&raw);
                        assert_eq!(arena.to_bytes(first), raw);
                        let independent = arena.from_bytes(&raw);
                        let mut changed = raw.clone();
                        changed[0] ^= 1;
                        let early = arena.from_bytes(&changed);
                        changed[0] ^= 1;
                        changed[size - 1] ^= 1;
                        let late = arena.from_bytes(&changed);
                        let (prefix, _) = arena.split(first, size as u64 - 1);
                        let ending = arena.from_bytes(&changed[size - 1..]);
                        let shared_late = arena.concat(prefix, ending);
                        let index = result.arenas.len();
                        for (label, second) in [("same_root", first), ("independent_equal", independent),
                            ("early_mismatch", early), ("late_mismatch", late), ("shared_late_mismatch", shared_late)] {
                            let bytes = arena.to_bytes(second);
                            assert_eq!(bytes.len(), size);
                            let expected = raw.iter().zip(&bytes).position(|(a,b)| a != b).unwrap_or(size) as u64;
                            assert_eq!(arena.lcp(first, second), expected);
                            assert_eq!(arena.equal(first, second), expected == size as u64);
                            result.fixtures.push(Fixture { arena: index, first, second, expected,
                                length: size as u64, name: format!("{size}/{label}") });
                        }
                        result.arenas.push(arena);
                    }
                    let mut arena = Arena::new();
                    let p = arena.from_bytes(b"ab");
                    let a = arena.repeat(p, 1u64 << 40);
                    let q = arena.from_bytes(b"abab");
                    let b = arena.repeat(q, 1u64 << 39);
                    assert_eq!(arena.len(a), 1u64 << 41);
                    assert_eq!(arena.lcp(a, b), 1u64 << 41);
                    assert!(arena.equal(a, b));
                    result.fixtures.push(Fixture { arena: 2, first: a, second: b,
                        expected: 1u64 << 41, length: 1u64 << 41, name: "repeat_2TiB".to_owned() });
                    let mut a = arena.from_bytes(b"ab");
                    let mut b = arena.from_bytes(b"ab");
                    for _ in 0..40 { a = arena.concat(a, a); b = arena.concat(b, b); }
                    assert_eq!(arena.len(a), 1u64 << 41);
                    assert_eq!(arena.lcp(a, b), 1u64 << 41);
                    assert!(arena.equal(a, b));
                    result.fixtures.push(Fixture { arena: 2, first: a, second: b,
                        expected: 1u64 << 41, length: 1u64 << 41, name: "self_concat_DAG_2TiB".to_owned() });
                    // Concrete default-profile collision guards the diagnostic's exact API.
                    let x = arena.from_bytes(&[0, 131]);
                    let y = arena.from_bytes(&[1, 0]);
                    assert_eq!(arena.hash(x), arena.hash(y));
                    assert_eq!(arena.lcp(x, y), 0);
                    assert!(!arena.equal(x, y));
                    result.arenas.push(arena);
                    result
                }
                pub fn len(&self) -> usize { self.fixtures.len() }
                pub fn name(&self, index: usize) -> &str { &self.fixtures[index].name }
                pub fn telemetry(&self, slot: &str) {
                    for fixture in &self.fixtures {
                        let arena = &self.arenas[fixture.arena];
                        let mut id = fixture.first.unwrap();
                        let root_address = arena.node(id) as *const NodeInner as usize;
                        let leaf_address = loop {
                            match arena.node(id) {
                                NodeInner::Leaf { data, .. } => break data.as_ptr() as usize,
                                NodeInner::Internal { left, .. } => id = *left,
                                NodeInner::Repeat { child, .. } => id = *child,
                            }
                        };
                        println!("fixture\t{slot}\t{}\t{}\t{}\t{}\t{}\t{}", fixture.name,
                            arena.node_count(), fixture.length, fixture.expected, root_address, leaf_address);
                    }
                }
                // Wrapper and operation dispatch are identical across slots. The operation
                // branch is outside the timed loop. Tiny cases include loop/wrapper overhead.
                #[inline(never)]
                pub fn measure(&self, index: usize, operation: Operation, iterations: usize) -> Row {
                    let fixture = &self.fixtures[index];
                    let arena = &self.arenas[fixture.arena];
                    let expected = match operation { Operation::Lcp => fixture.expected,
                        Operation::Equal => (fixture.expected == fixture.length) as u64 };
                    let mut answer = 0u64;
                    let processor_before = platform::processor();
                    let started = Instant::now();
                    match operation {
                        Operation::Lcp => for _ in 0..iterations {
                            answer = answer.wrapping_add(black_box(black_box(arena).lcp(black_box(fixture.first), black_box(fixture.second))));
                        },
                        Operation::Equal => for _ in 0..iterations {
                            answer = answer.wrapping_add(black_box(black_box(arena).equal(black_box(fixture.first), black_box(fixture.second))) as u64);
                        },
                    }
                    let elapsed_ns = started.elapsed().as_nanos();
                    let processor_after = platform::processor();
                    assert_eq!(answer, expected.wrapping_mul(iterations as u64));
                    Row { elapsed_ns, iterations, answer, processor_before, processor_after }
                }
            }
        }
    }
}
implementation!(a, slot_a);
implementation!(b, slot_b);

fn next(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut value = *seed;
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D049BB133111EB);
    value ^ (value >> 31)
}
fn main() {
    let arguments: Vec<String> = std::env::args().collect();
    let seed: u64 = arguments[1].parse().unwrap();
    let mask: usize = arguments[2].parse().unwrap();
    let block_ms: u64 = arguments.get(3).map(|x| x.parse().unwrap()).unwrap_or(10);
    let warm_ms: u64 = arguments.get(4).map(|x| x.parse().unwrap()).unwrap_or(10);
    let (initial, actual) = platform::affinity(mask);
    println!("meta\taffinity\t{initial}\t{actual}");
    println!("meta\tseed\t{seed}\t{block_ms}\t{warm_ms}");
    let (a, b) = if seed % 2 == 0 { let b = b::Suite::new(); let a = a::Suite::new(); (a,b) }
        else { (a::Suite::new(), b::Suite::new()) };
    assert_eq!(a.len(), b.len());
    a.telemetry("a"); b.telemetry("b");
    let mut order: Vec<usize> = (0..a.len()).collect();
    let mut random = seed;
    for i in (1..order.len()).rev() { order.swap(i, next(&mut random) as usize % (i + 1)); }
    for (case_order, index) in order.into_iter().enumerate() {
        assert_eq!(a.name(index), b.name(index));
        let operations = if (seed as usize + index) % 2 == 0 { [Operation::Lcp, Operation::Equal] }
            else { [Operation::Equal, Operation::Lcp] };
        for (api_order, operation) in operations.into_iter().enumerate() {
            // Independent of API position: each API gets ABBA and BAAB in
            // exactly three of the six processes, for every fixture.
            let first_a = (seed as usize + match operation { Operation::Lcp => 0, Operation::Equal => 1 }) % 2 == 0;
            // One shared calibration count; neither arm gets a different amount of work.
            let mut iterations = 16;
            let mut elapsed;
            loop {
                let (left, right) = if first_a {
                    (a.measure(index, operation, iterations), b.measure(index, operation, iterations))
                } else {
                    let right = b.measure(index, operation, iterations);
                    let left = a.measure(index, operation, iterations);
                    (left, right)
                };
                elapsed = left.elapsed_ns.max(right.elapsed_ns);
                if elapsed >= 1_000_000 { break; }
                iterations *= 2;
            }
            iterations = ((iterations as u128 * block_ms as u128 * 1_000_000) / elapsed).max(1) as usize;
            // Sustained paired warmup replaces the original three-operation warmup.
            let warm = Instant::now();
            let mut warm_a_first = first_a;
            while warm.elapsed() < Duration::from_millis(warm_ms) {
                if warm_a_first {
                    black_box(a.measure(index, operation, iterations));
                    black_box(b.measure(index, operation, iterations));
                } else {
                    black_box(b.measure(index, operation, iterations));
                    black_box(a.measure(index, operation, iterations));
                }
                warm_a_first = !warm_a_first;
            }
            let sequence = if first_a { ["a", "b", "b", "a"] }
                else { ["b", "a", "a", "b"] };
            for (position, slot) in sequence.into_iter().enumerate() {
                let row = if slot == "a" { a.measure(index, operation, iterations) }
                    else { b.measure(index, operation, iterations) };
                println!("time\t{}\t{}\t{case_order}\t{api_order}\t{position}\t{slot}\t{}\t{}\t{}\t{}\t{}",
                    a.name(index), operation.name(), row.iterations, row.elapsed_ns, row.answer,
                    row.processor_before, row.processor_after);
            }
        }
    }
}
