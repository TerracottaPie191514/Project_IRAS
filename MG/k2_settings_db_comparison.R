
prebuilt = import_biom("orig.biom")
custom_conf = import_biom("pb_conf.biom")
custom_no_conf = import_biom("pb_no_conf.biom")

prebuilt
custom_conf
custom_no_conf

# set rank names
colnames(tax_table(prebuilt)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(custom_conf)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(custom_no_conf)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# filtering 
prebuilt <- subset_taxa(prebuilt, Domain!="k__Archaea")
prebuilt <- subset_taxa(prebuilt, Domain!="k__Viruses")
prebuilt <- subset_taxa(prebuilt, Domain!="k__Eukaryota")

custom_conf <- subset_taxa(custom_conf, Domain!="k__Archaea")
custom_conf <- subset_taxa(custom_conf, Domain!="k__Viruses")
custom_conf <- subset_taxa(custom_conf, Domain!="k__Eukaryota")

custom_no_conf <- subset_taxa(custom_no_conf, Domain!="k__Archaea")
custom_no_conf <- subset_taxa(custom_no_conf, Domain!="k__Viruses")
custom_no_conf <- subset_taxa(custom_no_conf, Domain!="k__Eukaryota")

prebuilt
custom_conf
custom_no_conf
