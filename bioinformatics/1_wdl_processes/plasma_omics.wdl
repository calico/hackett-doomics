workflow biom {
    
    File path_do_config = "<<PATH_TO_do_config.json>>"
    Boolean overwrite = true
    
    call globals_parse_json {
      input: path_do_config=path_do_config
    }
    
    call featurization_parse_json {
      input: path_do_config=path_do_config
    }
    
    call printStringArray {
      input: input_string_array = featurization_parse_json.featurization_ids
    }
    
    scatter (featurization_id in featurization_parse_json.featurization_ids) {
      call featurization {
        input: process_feature_id = featurization_id,
          path_do_config = path_do_config,
          repo_path = globals_parse_json.repo_path,
          json_path_table = featurization_parse_json.featurization_table,
          overwrite = overwrite
      }
    }
}

task globals_parse_json {

    String path_do_config
    
    command <<<
        python <<CODE
        import json
        
        do_json_featurization = json.loads(open("${path_do_config}").read())["globals"]
        
        globals_dict = dict()
        for global_key in do_json_featurization.keys():
            global_key = str(global_key)
            global_value = do_json_featurization[global_key]
            globals_dict[global_key] = global_value
            print(global_key + ' = ' + global_value)
            
            with open(global_key + ".tsv", 'w') as out_h:
              out_h.write(global_value)
        CODE
    >>>    
    
    output {
      String project_path = read_string("project_path.tsv")
      String repo_path = read_string("repo_path.tsv")
      String do_drive_token = read_string("do_drive_token.tsv")
      String shock_token = read_string("shock_token.tsv")
    }
}

task featurization_parse_json {
    # takes a configuration file and determines the relevant featurizations that need to be generated
    
    String path_do_config
    
    command <<<
        python <<CODE
        import json
        
        do_json_featurization = json.loads(open("${path_do_config}").read())["featurization"]
        
        feature_dict = dict()
        for feature_type in do_json_featurization.keys():
            feature_type = str(feature_type)
            
            dataset_type_dict = dict()
            for dataset_type in do_json_featurization[feature_type].keys():
                dataset_type = str(dataset_type)
            
                dataset_type_keys = list(do_json_featurization[feature_type][dataset_type].keys())
                if len(dataset_type_keys) != 0:
                    dataset_type_dict[dataset_type] = dataset_type_keys
            
            if len(dataset_type_dict) != 0:
                feature_dict[feature_type] = dataset_type_dict
        
        print(feature_dict)
        
        # write featurization types as a tsv
        present_featurizations = list(feature_dict.keys())
        with open("featurizations.tsv", 'w') as out_h:
            out_h.write('\n'.join(present_featurizations))
            
        # write featurization dict as a tsv
        header = ["run_id", "featurization", "datatype", "runspecs"]

        with open("featurization_table.tsv", 'w') as out_h:
            out_h.write('\t'.join(header) + '\n')
    
            for f in feature_dict.keys():
                f = str(f)
                for d in feature_dict[f]:
                    d = str(d)
                    for r in feature_dict[f][d]:
                        r = str(r)
                        out_h.write('\t'.join([f + '.' + d + '.' + r, f, d, r]) + '\n')
        
        # create a map
        with open("featurization_map.tsv", 'w') as out_h:
            
            for f in feature_dict.keys():
                f = str(f)
                for d in feature_dict[f]:
                    d = str(d)
                    for r in feature_dict[f][d]:
                        r = str(r)
                        out_h.write('\t'.join([f + '.' + d + '.' + r, f]) + '\n')
        
        # unique features
            
        feature_ids = []
        for f in feature_dict.keys():
            f = str(f)
            for d in feature_dict[f]:
                d = str(d)
                for r in feature_dict[f][d]:
                    r = str(r)
                    feature_ids.append(f + '.' + d + '.' + r)
        
        with open("featurization_ids.tsv", 'w') as out_h:
            out_h.write('\n'.join(feature_ids))
        
        print("Featurization Data Types Parsed")
        CODE
    >>>
    
    output {
       Array[String] featurization_ids = read_lines("featurization_ids.tsv")
       Array[String] featurizations = read_lines("featurizations.tsv")
       Map[String, String] featurization_map = read_map("featurization_map.tsv")
       File featurization_table = "featurization_table.tsv"
    }
}

task featurization {
  
  String process_feature_id
  String path_do_config
  String repo_path
  File json_path_table
  Boolean overwrite
  
  String Rmd_rel_path = "/bioinformatics/utils/call_Rmd.R"
  String call_Rmd_path = repo_path + Rmd_rel_path
  
  command <<<
      Rscript ${call_Rmd_path} process_feature_id=${process_feature_id} json_path_table=${json_path_table} path_do_config=${path_do_config} overwrite=${overwrite}
  >>>   
    
  output { File out = "outputA.ext" }
}

task printStringArray {
  Array[String] input_string_array
  command {
    echo -e "input_string_array:\t${sep=',' input_string_array}"
  }
  output {
    String response = stdout()
  }
}