//! Golden outputs generated from the released, checksum-verified PyPI wheel.
//! All integers in the fixture are decimal strings, never JSON floating point.

use serde_json::Value;

fn fixture() -> Value {
    serde_json::from_str(include_str!("fixtures/python-0.2.2.json")).unwrap()
}

fn text<'a>(value: &'a Value, key: &str) -> &'a str {
    value[key].as_str().unwrap()
}

fn rows<'a>(value: &'a Value, key: &str) -> &'a [Value] {
    value[key].as_array().unwrap()
}

fn hex(data: &str) -> Vec<u8> {
    assert_eq!(data.len() % 2, 0);
    (0..data.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&data[i..i + 2], 16).unwrap())
        .collect()
}

/// Adapt only ownership and integer representation, keeping the test sequence
/// identical for both APIs and the Python oracle.
trait Engine {
    type Number;
    type Node: Copy;
    fn number(value: &str) -> Option<Self::Number>;
    fn new(prime: &Self::Number, base: &Self::Number) -> Self;
    fn build_bytes(&mut self, bytes: &[u8]) -> Self::Node;
    fn to_bytes(&self, node: Self::Node) -> Vec<u8>;
    fn hash(&self, node: Self::Node) -> String;
    fn len(&self, node: Self::Node) -> String;
    fn direct(&self, bytes: &[u8]) -> String;
    fn split(&mut self, node: Self::Node, at: &Self::Number) -> (Self::Node, Self::Node);
    fn concat(&mut self, left: Self::Node, right: Self::Node) -> Self::Node;
    fn repeat(&mut self, node: Self::Node, count: &Self::Number) -> Self::Node;
    fn substr(&mut self, node: Self::Node, start: &Self::Number, length: &Self::Number) -> String;
    fn nodes(&self) -> usize;
    fn validate(&self, node: Self::Node);
    fn power(&mut self, exponent: &Self::Number) -> String;
    fn phi(&self, count: &Self::Number, alpha: &Self::Number) -> String;
    fn modulo(&self, input: &str) -> Option<String>;
    fn overlap(&mut self, pattern: &[u8], length: &Self::Number) -> String;
}

struct Narrow(hashrope::Arena);

impl Engine for Narrow {
    type Number = u64;
    type Node = hashrope::Node;
    fn number(value: &str) -> Option<u64> {
        value.parse().ok()
    }
    fn new(prime: &u64, base: &u64) -> Self {
        Self(hashrope::Arena::with_hash(*prime, *base))
    }
    fn build_bytes(&mut self, bytes: &[u8]) -> Self::Node {
        self.0.from_bytes(bytes)
    }
    fn to_bytes(&self, node: Self::Node) -> Vec<u8> {
        self.0.to_bytes(node)
    }
    fn hash(&self, node: Self::Node) -> String {
        self.0.hash(node).to_string()
    }
    fn len(&self, node: Self::Node) -> String {
        self.0.len(node).to_string()
    }
    fn direct(&self, bytes: &[u8]) -> String {
        self.0.hash_bytes(bytes).to_string()
    }
    fn split(&mut self, node: Self::Node, at: &u64) -> (Self::Node, Self::Node) {
        self.0.split(node, *at)
    }
    fn concat(&mut self, left: Self::Node, right: Self::Node) -> Self::Node {
        self.0.concat(left, right)
    }
    fn repeat(&mut self, node: Self::Node, count: &u64) -> Self::Node {
        self.0.repeat(node, *count)
    }
    fn substr(&mut self, node: Self::Node, start: &u64, length: &u64) -> String {
        self.0.substr_hash(node, *start, *length).to_string()
    }
    fn nodes(&self) -> usize {
        self.0.node_count()
    }
    fn validate(&self, node: Self::Node) {
        self.0.validate(node)
    }
    fn power(&mut self, exponent: &u64) -> String {
        self.0.hasher_mut().power(*exponent).to_string()
    }
    fn phi(&self, count: &u64, alpha: &u64) -> String {
        hashrope::phi(*count, *alpha, self.0.hasher().prime()).to_string()
    }
    fn modulo(&self, input: &str) -> Option<String> {
        input
            .parse::<u128>()
            .ok()
            .map(|a| hashrope::mersenne_mod(a, self.0.hasher().prime()).to_string())
    }
    fn overlap(&mut self, pattern: &[u8], length: &u64) -> String {
        let h = self.0.hasher_mut();
        let remainder = (*length % pattern.len() as u64) as usize;
        h.hash_overlap(
            h.hash(pattern),
            pattern.len() as u64,
            *length,
            h.hash(&pattern[..remainder]),
        )
        .to_string()
    }
}

fn assert_node<E: Engine>(engine: &E, node: E::Node, expected: &Value) {
    assert_eq!(engine.hash(node), text(expected, "hash"));
    assert_eq!(engine.len(node), text(expected, "length"));
    engine.validate(node);
}

fn assert_ranges<E: Engine>(engine: &mut E, node: E::Node, ranges: &[Value]) {
    for range in ranges {
        let start = E::number(text(range, "start")).unwrap();
        let length = E::number(text(range, "length")).unwrap();
        let nodes_before = engine.nodes();
        assert_eq!(engine.substr(node, &start, &length), text(range, "hash"));
        assert_eq!(
            engine.nodes(),
            nodes_before,
            "substring hashing allocated rope nodes"
        );
    }
}

fn run_profiles<E: Engine>() -> usize {
    let vectors = fixture();
    assert_eq!(text(&vectors["reference"], "version"), "0.2.2");
    let mut tested_profiles = 0;
    for profile in rows(&vectors, "profiles") {
        let Some(prime) = E::number(text(profile, "prime")) else {
            continue;
        };
        let base = E::number(text(profile, "base")).unwrap();
        let mut engine = E::new(&prime, &base);
        for row in rows(profile, "bytes") {
            let bytes = hex(text(row, "bytes_hex"));
            let node = engine.build_bytes(&bytes);
            assert_node(&engine, node, row);
            assert_eq!(engine.direct(&bytes), text(row, "hash"));
            assert_eq!(engine.to_bytes(node), bytes);
            for split in rows(row, "splits") {
                let cut = text(split, "at");
                let at = E::number(cut).unwrap();
                let (left, right) = engine.split(node, &at);
                assert_node(&engine, left, &split["left"]);
                assert_node(&engine, right, &split["right"]);
                let index: usize = cut.parse().unwrap();
                assert_eq!(engine.to_bytes(left), bytes[..index]);
                assert_eq!(engine.to_bytes(right), bytes[index..]);
                let joined = engine.concat(left, right);
                assert_node(&engine, joined, row);
                // Persistence: edits must not alter a previously built root.
                assert_eq!(engine.to_bytes(node), bytes);
            }
            assert_ranges(&mut engine, node, rows(row, "substrings"));
        }
        for row in rows(profile, "repeats") {
            let Some(count) = E::number(text(row, "count")) else {
                continue;
            };
            if E::number(text(row, "length")).is_none() {
                continue;
            }
            let child = engine.build_bytes(&hex(text(row, "pattern_hex")));
            let node = engine.repeat(child, &count);
            assert_node(&engine, node, row);
            assert_ranges(&mut engine, node, rows(row, "substrings"));
            if text(row, "count").parse::<usize>().is_ok_and(|n| n < 100) {
                let count: usize = text(row, "count").parse().unwrap();
                assert_eq!(
                    engine.to_bytes(node),
                    hex(text(row, "pattern_hex")).repeat(count)
                );
            }
        }
        for row in rows(profile, "powers") {
            if let Some(exponent) = E::number(text(row, "exponent")) {
                assert_eq!(engine.power(&exponent), text(row, "value"));
                assert_eq!(engine.power(&exponent), text(row, "value"), "cached power");
            }
        }
        for row in rows(profile, "phis") {
            if let Some(count) = E::number(text(row, "count")) {
                let alpha = E::number(text(row, "alpha")).unwrap();
                assert_eq!(engine.phi(&count, &alpha), text(row, "value"));
            }
        }
        for row in rows(profile, "moduli") {
            if let Some(value) = engine.modulo(text(row, "input")) {
                assert_eq!(value, text(row, "value"));
            }
        }
        for row in rows(profile, "overlaps") {
            if let Some(length) = E::number(text(row, "length")) {
                assert_eq!(
                    engine.overlap(&hex(text(row, "pattern_hex")), &length),
                    text(row, "hash")
                );
            }
        }
        let edits = rows(profile, "edits");
        let original = engine.build_bytes(&hex(text(&edits[0], "bytes_hex")));
        let mut node = original;
        assert_node(&engine, node, &edits[0]);
        for row in &edits[1..] {
            let start = E::number(text(row, "start")).unwrap();
            let removed = E::number(text(row, "removed")).unwrap();
            let (left, rest) = engine.split(node, &start);
            let (_, right) = engine.split(rest, &removed);
            let insert = engine.build_bytes(&hex(text(row, "insert_hex")));
            let prefix = engine.concat(left, insert);
            node = engine.concat(prefix, right);
            assert_node(&engine, node, row);
            assert_eq!(engine.to_bytes(node), hex(text(row, "bytes_hex")));
            assert_node(&engine, original, &edits[0]);
        }
        tested_profiles += 1;
    }
    tested_profiles
}

#[test]
fn published_python_valid_byte_contract_u64() {
    assert_eq!(run_profiles::<Narrow>(), 5);
}

#[test]
fn published_python_streams_legacy_and_bounded() {
    for profile in rows(&fixture(), "profiles") {
        let Ok(prime) = text(profile, "prime").parse::<u64>() else {
            continue;
        };
        let base = text(profile, "base").parse().unwrap();
        for bounded in [false, true] {
            let mut stream = if bounded {
                hashrope::SlidingWindow::new_bounded(8, 4, prime, base)
            } else {
                hashrope::SlidingWindow::new(8, 4, prime, base)
            };
            for op in rows(profile, "stream") {
                match text(op, "kind") {
                    "bytes" => stream.append_bytes(&hex(text(op, "bytes_hex"))),
                    "copy" => stream.append_copy(
                        text(op, "offset").parse().unwrap(),
                        text(op, "length").parse().unwrap(),
                    ),
                    _ => unreachable!(),
                }
                assert_eq!(stream.current_hash().to_string(), text(op, "hash"));
                assert_eq!(stream.final_hash().to_string(), text(op, "hash"));
                assert_eq!(stream.pos().to_string(), text(op, "position"));
                assert_eq!(stream.window_len().to_string(), text(op, "window_length"));
            }
        }
    }
}

#[test]
fn lazy_default_profile_preserves_sentinel_and_materializes_python_hashes() {
    let vectors = fixture();
    let profile = &vectors["profiles"][3];
    let mut arena = hashrope::Arena::new_lazy();
    for row in rows(profile, "bytes") {
        let node = arena.from_bytes(&hex(text(row, "bytes_hex")));
        if let Some(id) = node {
            assert_eq!(arena.hash(node), hashrope::LAZY_SENTINEL);
            arena.ensure_hash(id);
        }
        assert_eq!(arena.hash(node).to_string(), text(row, "hash"));
        for sub in rows(row, "substrings") {
            assert_eq!(
                arena
                    .substr_hash(
                        node,
                        text(sub, "start").parse().unwrap(),
                        text(sub, "length").parse().unwrap()
                    )
                    .to_string(),
                text(sub, "hash")
            );
        }
    }
}

#[cfg(feature = "bigint")]
mod wide {
    use super::*;
    use hashrope::bigint::{self, BigUint};
    struct Wide(bigint::Arena);

    impl Engine for Wide {
        type Number = BigUint;
        type Node = bigint::Node;
        fn number(value: &str) -> Option<BigUint> {
            value.parse().ok()
        }
        fn new(prime: &BigUint, base: &BigUint) -> Self {
            Self(bigint::Arena::with_hash(prime.clone(), base.clone()))
        }
        fn build_bytes(&mut self, bytes: &[u8]) -> Self::Node {
            self.0.from_bytes(bytes)
        }
        fn to_bytes(&self, node: Self::Node) -> Vec<u8> {
            self.0.to_bytes(node)
        }
        fn hash(&self, node: Self::Node) -> String {
            self.0.hash(node).to_string()
        }
        fn len(&self, node: Self::Node) -> String {
            self.0.len(node).to_string()
        }
        fn direct(&self, bytes: &[u8]) -> String {
            self.0.hash_bytes(bytes).to_string()
        }
        fn split(&mut self, node: Self::Node, at: &BigUint) -> (Self::Node, Self::Node) {
            self.0.split(node, at)
        }
        fn concat(&mut self, left: Self::Node, right: Self::Node) -> Self::Node {
            self.0.concat(left, right)
        }
        fn repeat(&mut self, node: Self::Node, count: &BigUint) -> Self::Node {
            self.0.repeat(node, count)
        }
        fn substr(&mut self, node: Self::Node, start: &BigUint, length: &BigUint) -> String {
            self.0.substr_hash(node, start, length).to_string()
        }
        fn nodes(&self) -> usize {
            self.0.node_count()
        }
        fn validate(&self, node: Self::Node) {
            self.0.validate(node)
        }
        fn power(&mut self, exponent: &BigUint) -> String {
            self.0.hasher_mut().power(exponent).to_string()
        }
        fn phi(&self, count: &BigUint, alpha: &BigUint) -> String {
            bigint::phi(count, alpha, self.0.hasher().prime()).to_string()
        }
        fn modulo(&self, input: &str) -> Option<String> {
            Some(bigint::mersenne_mod(&input.parse().unwrap(), self.0.hasher().prime()).to_string())
        }
        fn overlap(&mut self, pattern: &[u8], length: &BigUint) -> String {
            let h = self.0.hasher_mut();
            let size = BigUint::from(pattern.len());
            let remainder: usize = (length % &size).try_into().unwrap();
            h.hash_overlap(
                &h.hash(pattern),
                &size,
                length,
                &h.hash(&pattern[..remainder]),
            )
            .to_string()
        }
    }

    #[test]
    fn published_python_arbitrary_precision_contract() {
        assert_eq!(run_profiles::<Wide>(), 8);
    }

    #[test]
    fn published_python_arbitrary_precision_streams() {
        for profile in rows(&fixture(), "profiles") {
            let prime: BigUint = text(profile, "prime").parse().unwrap();
            let base: BigUint = text(profile, "base").parse().unwrap();
            for bounded in [false, true] {
                let mut stream = if bounded {
                    bigint::SlidingWindow::new_bounded(
                        8u8.into(),
                        4u8.into(),
                        prime.clone(),
                        base.clone(),
                    )
                } else {
                    bigint::SlidingWindow::new(8u8.into(), 4u8.into(), prime.clone(), base.clone())
                };
                for op in rows(profile, "stream") {
                    match text(op, "kind") {
                        "bytes" => stream.append_bytes(&hex(text(op, "bytes_hex"))),
                        "copy" => stream.append_copy(
                            &text(op, "offset").parse().unwrap(),
                            &text(op, "length").parse().unwrap(),
                        ),
                        _ => unreachable!(),
                    }
                    assert_eq!(stream.current_hash().to_string(), text(op, "hash"));
                    assert_eq!(stream.final_hash().to_string(), text(op, "hash"));
                    assert_eq!(stream.pos().to_string(), text(op, "position"));
                    assert_eq!(stream.window_len().to_string(), text(op, "window_length"));
                }
            }
        }
    }
}

#[test]
fn published_pollard_jsonl_remains_byte_and_hash_stable() {
    use sha2::{Digest, Sha256};
    let log = include_bytes!("fixtures/pollard-1.6.0-hashrope-0.2.2.jsonl");
    assert_eq!(log.len(), 785);
    assert_eq!(
        format!("{:x}", Sha256::digest(log)),
        "45174ce9b17341c5bd9efb800a403308d148ce799b60862896341c248a37e623"
    );
    let expected_hash = 2_269_240_639_182_000_324;
    let expected_lines = [
        1_416_890_852_347_916_452,
        937_132_792_245_618_969,
        2_174_710_591_043_672_813,
    ];
    let mut arena = hashrope::Arena::new();
    let rebuilt = arena.from_bytes(log);
    let mut incremental = None;
    for (line, expected_line) in log.split_inclusive(|&b| b == b'\n').zip(expected_lines) {
        let part = arena.from_bytes(line);
        assert_eq!(arena.hash(part), expected_line);
        incremental = arena.concat(incremental, part);
    }
    assert_eq!(arena.hash_bytes(log), expected_hash);
    assert_eq!(arena.hash(rebuilt), expected_hash);
    assert_eq!(arena.hash(incremental), expected_hash);
    assert_eq!(arena.to_bytes(incremental), log);
    arena.validate(incremental);
    // A framework may stream transport chunks through arbitrary UTF-8/JSON
    // boundaries. Neither JSON normalization nor text decoding is performed.
    let mut chunks = None;
    let mut window = hashrope::SlidingWindow::new_bounded(32, 8, hashrope::MERSENNE_61, 131);
    for bytes in log.chunks(7) {
        let part = arena.from_bytes(bytes);
        chunks = arena.concat(chunks, part);
        window.append_bytes(bytes);
    }
    assert_eq!(arena.hash(chunks), expected_hash);
    assert_eq!(window.final_hash(), expected_hash);
}
