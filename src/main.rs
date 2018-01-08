mod utils;
mod configmaker;
mod variables;
mod parameter;
mod fcalculator;
mod meshlist;
mod observer;
mod mdsystem;
use mdsystem::MDSystem;

fn main() {
    let mut mdsystem = MDSystem::new();
    mdsystem.run();
}

