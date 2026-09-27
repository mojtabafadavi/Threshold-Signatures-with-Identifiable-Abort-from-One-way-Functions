use rand::thread_rng;
use std::time::Instant;
use tlots::shamir;

fn main() {
    println!("Empirical Dist cost: share_naive (O(N*T)) vs share_ntt (O(N log N))");
    println!("{:>8} {:>8} {:>16} {:>16} {:>10}", "N", "T", "naive (ms)", "ntt (ms)", "speedup");
    let mut rng = thread_rng();
    let secret = 123456789u64 % shamir::P;

    // Naive is O(N*T); cap it here so the comparison finishes in reasonable time.
    for &n in &[100usize, 500, 1000, 2000, 4000, 8000] {
        let t = n / 2;

        let start = Instant::now();
        let _naive = shamir::share_naive(secret, t, n, &mut rng);
        let naive_time = start.elapsed().as_secs_f64() * 1000.0;

        let start = Instant::now();
        let _ntt = shamir::share_ntt(secret, t, n, &mut rng);
        let ntt_time = start.elapsed().as_secs_f64() * 1000.0;

        println!(
            "{:>8} {:>8} {:>16.3} {:>16.3} {:>9.1}x",
            n, t, naive_time, ntt_time, naive_time / ntt_time
        );
    }

    println!("\nNTT alone, scaled far beyond where naive is even attemptable:");
    println!("{:>8} {:>8} {:>16}", "N", "T", "ntt (ms)");
    for &n in &[16000usize, 32000, 65536, 131072, 262144, 524288] {
        let t = n / 2;
        let start = Instant::now();
        let _ntt = shamir::share_ntt(secret, t, n, &mut rng);
        let ntt_time = start.elapsed().as_secs_f64() * 1000.0;
        println!("{:>8} {:>8} {:>16.3}", n, t, ntt_time);
    }
}
