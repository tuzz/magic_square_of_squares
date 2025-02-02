use clap::Parser;

#[derive(Parser)]
pub struct Args {
    #[arg(long, help = "The magic sum to start checking from. Must be congruent to 3 modulo 72.", default_value = "3")]
    pub from: u64,

    #[arg(long, help = "The amount to increase the magic sum by each time. Must be divisible by 72.", default_value = "72")]
    pub stride: u64,

    #[arg(long, help = "The number of checks to perform in parallel using operating system threads.", default_value_t = num_threads())]
    pub threads: u64
}

impl Args {
    pub fn parse_command_line() -> Self {
        let args = Self::parse();

        if args.from % 72 != 3 {
            eprintln!("Error: Must start from a number that is 3 mod 72 such as 3, 75, 147, etc.");
            std::process::exit(1);
        }

        if args.stride == 0 || args.stride % 72 != 0 {
            eprintln!("Error: Must stride by a number that is divisible by 72 such as 72, 144, 216, etc.");
            std::process::exit(1);
        }

        if args.threads == 0 {
            eprintln!("Error: Must use at least one thread. Omit --threads to use the system default.");
            std::process::exit(1);
        }

        args
    }
}

fn num_threads() -> u64 {
    std::thread::available_parallelism().map(|n| n.into()).unwrap_or(1) as u64
}
