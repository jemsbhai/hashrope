//! New workload, not the historical EXP020 fingerprint-read benchmark.
//! 50 seeded ASCII churn edits, each followed by comparison to a prebuilt target.
//! Targets alternate equal / final-byte-different (50% each). Independent and
//! shared-history targets are reported separately; target construction is untimed.
use hashrope::{Arena, Node};
use ropey::Rope;
use std::{hint::black_box, time::Instant};

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
fn edits(n: usize, seed: u64) -> (Vec<u8>, Vec<Edit>, Vec<Vec<u8>>) {
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
        targets.push(target);
        operations.push(edit);
    }
    (initial, operations, targets)
}
fn edit_hashrope(a: &mut Arena, root: Node, edit: &Edit) -> Node {
    let (left, rest) = a.split(root, edit.delete as u64);
    let (_, right) = a.split(rest, edit.bytes.len() as u64);
    let root = a.concat(left, right);
    let (left, right) = a.split(root, edit.insert as u64);
    let added = a.from_bytes(&edit.bytes);
    let head = a.concat(left, added);
    a.concat(head, right)
}
fn edit_ropey(rope: &mut Rope, edit: &Edit) {
    rope.remove(edit.delete..edit.delete + edit.bytes.len());
    rope.insert(edit.insert, std::str::from_utf8(&edit.bytes).unwrap());
}
fn main() {
    let seed: u64 = std::env::args().nth(1).unwrap().parse().unwrap();
    let invocation: usize = std::env::args().nth(2).unwrap_or_else(|| "0".to_string()).parse().unwrap();
    for n in [10_000usize, 100_000, 1_000_000] {
        let (initial, operations, expected_targets) = edits(n, seed);
        for shared in [false, true] {
            let mut arms = ["hashrope_exact", "hashrope_candidate", "ropey_exact", "ropey_rehash_candidate"];
            arms.rotate_left((seed as usize + invocation) % 4);
            for arm in arms {
                for sample in 0..3 {
                    let elapsed;
                    let mut equal_count = 0;
                    let mut equal_mask = 0u64;
                    if arm.starts_with("hashrope") {
                        let mut arena = Arena::new();
                        let mut root = arena.from_bytes(&initial);
                        let mut independent = Arena::new();
                        let mut targets = Vec::new();
                        if shared {
                            let mut target_root = root;
                            for (i, edit) in operations.iter().enumerate() {
                                target_root = edit_hashrope(&mut arena, target_root, edit);
                                let target = if i % 2 == 1 {
                                    let (prefix, _) = arena.split(target_root, n as u64 - 1);
                                    let ending = arena.from_bytes(&expected_targets[i][n - 1..]);
                                    arena.concat(prefix, ending)
                                } else { target_root };
                                targets.push(target);
                            }
                        } else {
                            for target in &expected_targets { targets.push(independent.from_bytes(target)); }
                        }
                        // Warm representative cache entries with untimed edits.
                        let mut warm = root;
                        for edit in &operations { warm = edit_hashrope(&mut arena, warm, edit); }
                        black_box(warm);
                        let start = Instant::now();
                        for (i, edit) in operations.iter().enumerate() {
                            root = edit_hashrope(&mut arena, root, edit);
                            let equal = if arm == "hashrope_exact" {
                                if shared { arena.equal(root, targets[i]) }
                                else { arena.equal_across(root, &independent, targets[i]) }
                            } else if shared { arena.hash(root) == arena.hash(targets[i]) }
                            else { arena.hash(root) == independent.hash(targets[i]) };
                            equal_count += black_box(equal) as usize;
                            equal_mask |= (equal as u64) << i;
                        }
                        elapsed = start.elapsed();
                        assert_eq!(arena.to_bytes(root), {
                            // Last target is a mismatch; verify against deterministic edits directly.
                            let mut last = initial.clone();
                            for edit in &operations {
                                last.drain(edit.delete..edit.delete + edit.bytes.len());
                                last.splice(edit.insert..edit.insert, edit.bytes.iter().copied());
                            }
                            last
                        });
                    } else {
                        let mut rope = Rope::from_str(std::str::from_utf8(&initial).unwrap());
                        let mut targets = Vec::new();
                        if shared {
                            let mut target = rope.clone();
                            for (i, edit) in operations.iter().enumerate() {
                                edit_ropey(&mut target, edit);
                                let mut copy = target.clone();
                                if i % 2 == 1 {
                                    copy.remove(n - 1..n);
                                    copy.insert(n - 1, std::str::from_utf8(&expected_targets[i][n - 1..]).unwrap());
                                }
                                targets.push(copy);
                            }
                        } else {
                            for target in &expected_targets { targets.push(Rope::from_str(std::str::from_utf8(target).unwrap())); }
                        }
                        let hasher = Arena::new();
                        let target_hashes: Vec<u64> = expected_targets.iter().map(|t| hasher.hash_bytes(t)).collect();
                        let start = Instant::now();
                        for (i, edit) in operations.iter().enumerate() {
                            edit_ropey(&mut rope, edit);
                            let equal = if arm == "ropey_exact" { rope == targets[i] }
                            else { hasher.hash_bytes(rope.to_string().as_bytes()) == target_hashes[i] };
                            equal_count += black_box(equal) as usize;
                            equal_mask |= (equal as u64) << i;
                        }
                        elapsed = start.elapsed();
                    }
                    assert_eq!(equal_count, 25, "{arm}, n={n}, shared={shared}");
                    assert_eq!(equal_mask, (0..50).step_by(2).map(|i| 1u64 << i).sum::<u64>());
                    println!("{seed}\t{n}\t{shared}\t{arm}\t{sample}\t{}", elapsed.as_secs_f64() * 1000.0 / 50.0);
                }
            }
        }
    }
}
