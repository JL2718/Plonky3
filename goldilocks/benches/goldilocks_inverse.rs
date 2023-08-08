use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use p3_field::{AbstractField, Field, PrimeField32,inverse_u, try_inverse_u_loop, PrimeField64, try_inverse_u64,try_inverse_u64_loop};
use p3_goldilocks::Goldilocks;

type F = Goldilocks;

fn try_inverse(c: &mut Criterion) {
    c.bench_function("power_inverse", |b| {
        b.iter_batched(
            || F::from_canonical_u64(rand::random::<u64>()),
            |x| x.try_inverse(),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("power_inverse_2", |b| {
        b.iter_batched(
            || F::from_canonical_u64(rand::random::<u64>()),
            |x| x.exp_u64(F::ORDER_U64-2),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("gcd_loop_inverse", |b| {
        b.iter_batched(
            || F::from_canonical_u64(rand::random::<u64>()).as_canonical_u64(),
            |x| inverse_u::<u64, i64>(x, F::ORDER_U64),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("gcd_inverse", |b| {
        b.iter_batched(
            || F::from_canonical_u64(rand::random::<u64>()).as_canonical_u64(),
            |x| try_inverse_u_loop::<u64, i64>(x, F::ORDER_U64),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("gcd_inverse_u64", |b| {
        b.iter_batched(
            || F::from_canonical_u64(rand::random::<u64>()).as_canonical_u64(),
            |x| try_inverse_u64(x, F::ORDER_U64),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("gcd_inverse_u64_loop", |b| {
        b.iter_batched(
            || F::from_canonical_u64(rand::random::<u64>()).as_canonical_u64(),
            |x| try_inverse_u64_loop(x, F::ORDER_U64),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(goldilocks_arithmetics, try_inverse);
criterion_main!(goldilocks_arithmetics);
