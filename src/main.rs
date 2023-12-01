use std::fs;
use aho_corasick::AhoCorasick;
use csv::ReaderBuilder;
use std::fs::File;
use std::collections::HashMap;
use std::io;
use plotters::prelude::*;


fn main() -> std::io::Result<()> {
    let (uniprot_vec, protein_seq_vec) = fasta_read("UniProt_Human.fasta")?;
    let peptide_vec = psv_read("psm.tsv");
    let _match_hash = aho_search(protein_seq_vec, peptide_vec?, uniprot_vec);
    drawings();
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

fn drawings() -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("plotters-doc-data/0.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-1f32..1f32, -0.1f32..1f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            (-50..=50).map(|x| x as f32 / 50.0).map(|x| (x, x * x)),
            &RED,
        ))?
        .label("y = x^2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}