use std::fs;
use aho_corasick::{AhoCorasick, PatternID};
use csv::ReaderBuilder;
use std::fs::File;
use std::collections::HashMap;
use std::io;
use std::io::Write;
use bokeh::{params::KERNEL9_PARAM_SET, Blur};
use image::{io::Reader as ImageReader, ImageError};


fn main() -> std::io::Result<()> {
    let (uniprot_vec, protein_seq_vec) = fasta_read("UniProt_Human.fasta")?;
    let peptide_vec = psv_read("psm.tsv");
    let match_hash = aho_search(protein_seq_vec, peptide_vec?, uniprot_vec);
    drawing();
    Ok(())
}

fn fasta_read(file_name: &str) -> std::io::Result<(Vec<String>, Vec<String>)> {
    // read the uniprot human fasta file
    let file_content = fs::read_to_string(file_name)?;
    // split vector contains the individual proteins split with \n>
    let split_vec: Vec<&str> = file_content.split("\n>").collect();
    let uniprot_vec: Vec<String> = split_vec.iter().filter_map(|&x| x.split("|").nth(1).map(|s| s.to_string())).collect();
    let _protein_seq_vec: Vec<String> = split_vec.iter().filter_map(|&x| {let x: Vec<&str> = x.split('\n').collect();if x.len() > 1 {Some(x[1..].join(""))} else {None}}).collect();
    Ok((uniprot_vec, _protein_seq_vec))
}

fn psv_read(file_name: &str)  -> std::io::Result<Vec<String>> {
    let file = File::open(file_name)?;
    let mut tsv_rdr: csv::Reader<_> = ReaderBuilder::new().delimiter(b'\t').from_reader(file);
    let mut peptides_vec:Vec<String> = Vec::new();
    for result in tsv_rdr.records() {
        let record = result?;
                if let Some(value) = record.get(2) {
            peptides_vec.push(value.to_string());
        }
    }

    Ok(peptides_vec)
    
}


fn aho_search(protein_seq_vec: Vec<String>, peptides_vec: Vec<String>, uniprot_vec: Vec<String>) -> io::Result<HashMap<String, Vec<String>>> {
    let search_seq: String = protein_seq_vec.join("|"); 
    let automaton = AhoCorasick::new(&peptides_vec).unwrap();

    let uniprot_id:Vec<String> = uniprot_vec.iter().zip(protein_seq_vec.iter())
        .flat_map(|(item, seq)| {
            let mut repeated_items = std::iter::repeat(item.clone())
                .take(seq.len())
                .collect::<Vec<_>>();
            repeated_items.push("|".to_string());
            repeated_items
        })
        .collect::<Vec<String>>();

    let mut match_dict: HashMap<String, Vec<String>> = HashMap::new();

    for mat in automaton.find_iter(&search_seq) {
        let uniprot_key: &String = &uniprot_id[mat.start()];
        let peptide_match: String = peptides_vec[mat.pattern()].clone();
        match_dict.entry(uniprot_key.clone()).or_insert_with(Vec::new).push(peptide_match);
    }

    Ok(match_dict)
}

fn drawing() -> std::io::Result<()> {
    let data = vec![1, 2, 3, 4, 5];
    let mut file = File::create("data.csv")?;
    
    writeln!(file, "x")?;
    for x in data {
        writeln!(file, "{}", x)?;
    }
    Ok(())
}