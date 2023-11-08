use std::io::prelude::*;
use std::fs::File;
use std::io::BufReader;

fn main() -> std::io::Result<()> {
    let fin = File::open("psm.tsv")?;
    let mut reader = BufReader::new(fin);

    let mut line = String::new();
    let len = reader.read_line(&mut line)?;
    println!("{len} bytes long");
    Ok(())
}
