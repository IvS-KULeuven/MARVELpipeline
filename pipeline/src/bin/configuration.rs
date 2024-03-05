

pub mod configuration {
    use std::fs::File;
    use std::env;
    pub fn parse_file () -> serde_yaml::Value {
        // Load and parse the param.yaml file to get the paths from which we
        // will load the files and to which we will save the output files.
        let path = env::current_dir().unwrap();
        let filename = path.join("../params.yaml");
        let filename = filename.as_path();

        let configuration: serde_yaml::Value = match File::open(filename) {
            Ok(file) => match serde_yaml::from_reader(file) {
                Ok(data) =>data,
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

        configuration
    }
}
