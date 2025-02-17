
use std::fs::File;
use std::env;


// Load and parse the param.yaml file to get the paths from which we
// will load the files and to which we will save the output files.

pub fn parse_file(config_path: &str) -> serde_yaml::Value {
    let path = env::current_dir().unwrap();
    let filename = path.join(config_path);
    let filename = filename.as_path();

    let config: serde_yaml::Value = match File::open(filename) {
        Ok(file) => match serde_yaml::from_reader(file) {
            Ok(data) => data,
            Err(err) => {
                eprintln!("There was an error parsing the YAML file{}", err);
                std::process::exit(1);
            }
        },
        Err(error) => {
            eprintln!("Error opening file {}: {}", filename.display(), error);
            std::process::exit(1);
        }
    };

    config
}
