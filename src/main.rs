mod utils;
mod configmaker;
mod variables;
mod parameter;
mod fcalculator;
mod meshlist;
mod observer;
mod mdsystem;
use mdsystem::MDSystem;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("Usage:");
        println!("$ {} input_dir", args[0]);
        panic!("directory path should be specified.");
    }
    let mut mdsystem = MDSystem::new(&args[1]);
    mdsystem.run();
}
