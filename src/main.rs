use std::fs;

fn main() -> std::io::Result<()> {
    // read the uniprot human fasta file
    let file_content = fs::read_to_string("UniProt_Human.fasta")?;
    // split vector contains the individual proteins split with \n>
    let split_vec: Vec<&str> = file_content.split("\n>").collect();
    println!("{}", &split_vec[31]);
    let uniprot_vec: Vec<&str> = split_vec.iter().filter_map(|&x| x.split("|").nth(1)).collect();
    println!("{}", &uniprot_vec[31]);
    let protein_seq_vec: Vec<&str> = split_vec.iter().filter_map(|&x| x.split("|").nth(1)).collect();
    Ok(())
}