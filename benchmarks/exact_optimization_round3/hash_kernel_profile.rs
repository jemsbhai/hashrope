//! Diagnostic only: compare mathematically equivalent default-profile kernels.
use std::{hint::black_box, time::Instant};
const P: u64 = (1u64 << 61) - 1;
#[inline]
fn reduce(a: u128, p: u64) -> u64 {
    let k = 64 - p.leading_zeros();
    let mut v = (a >> k) + (a & p as u128);
    v = (v >> k) + (v & p as u128);
    let r = v as u64;
    if r >= p { r - p } else { r }
}
#[inline(never)]
fn legacy(data: &[u8], p: u64, x: u64) -> u64 {
    let mut h = 0u64;
    for &b in data { h = reduce(h as u128 * x as u128 + b as u128 + 1, p); }
    h
}
#[inline(never)]
fn constant(data: &[u8], p: u64, x: u64) -> u64 {
    if p != P || x != 131 { return legacy(data, p, x); }
    let mut h = 0u64;
    for &b in data { h = reduce(h as u128 * 131 + b as u128 + 1, P); }
    h
}
#[inline(never)]
fn grouped<const N: usize, const DOT: bool>(data: &[u8], p: u64, x: u64) -> u64 {
    if p != P || x != 131 { return legacy(data, p, x); }
    let mut h = 0u64;
    let mut chunks = data.chunks_exact(N);
    for chunk in &mut chunks {
        let mut block = 0u64;
        if DOT {
            for (i, &b) in chunk.iter().enumerate() { block += (b as u64 + 1) * 131u64.pow((N - i - 1) as u32); }
        } else {
            for &b in chunk { block = block * 131 + b as u64 + 1; }
        }
        h = reduce(h as u128 * 131u64.pow(N as u32) as u128 + block as u128, P);
    }
    for &b in chunks.remainder() { h = reduce(h as u128 * 131 + b as u128 + 1, P); }
    h
}
fn main() {
    let cases: [(&str, fn(&[u8],u64,u64)->u64); 6] = [
        ("legacy",legacy), ("constant",constant), ("horner4",grouped::<4,false>),
        ("dot4",grouped::<4,true>), ("horner8",grouped::<8,false>), ("dot8",grouped::<8,true>)];
    for len in [0,1,7,8,16,31,64,255,256,511,512,4096,65536] {
        let data: Vec<u8> = (0..len).map(|i| (i*37+13) as u8).collect();
        let expected = legacy(&data,P,131);
        let iterations = (2_000_000 / (len+1)).max(100);
        for (name,f) in cases {
            assert_eq!(f(&data,P,131), expected);
            for _ in 0..100 { black_box(f(black_box(&data),black_box(P),black_box(131))); }
            for sample in 0..7 {
                let start = Instant::now();
                for _ in 0..iterations { black_box(f(black_box(&data),black_box(P),black_box(131))); }
                println!("{{\"kernel\":\"{name}\",\"length\":{len},\"sample\":{sample},\"ns\":{}}}",start.elapsed().as_nanos() as f64/iterations as f64);
            }
        }
    }
}
