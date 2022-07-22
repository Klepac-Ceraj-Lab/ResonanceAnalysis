const microbiome_predictors = [
    "Absiella_dolichum", "Acetobacter_cibinongensis", "Acidaminococcus_intestini", "Acidaminococcus_sp_CAG_542", "Acinetobacter_bereziniae",
    "Acinetobacter_johnsonii", "Acinetobacter_pittii", "Acinetobacter_soli", "Acinetobacter_ursingii", "Actinomyces_europaeus",
    "Actinomyces_graevenitzii", "Actinomyces_johnsonii", "Actinomyces_odontolyticus", "Actinomyces_oris", "Actinomyces_radingae",
    "Actinomyces_sp_HMSC035G02", "Actinomyces_sp_HPA0247", "Actinomyces_sp_ICM47", "Actinomyces_sp_S6_Spd3", "Actinomyces_sp_oral_taxon_181",
    "Actinomyces_turicensis", "Actinomyces_urogenitalis", "Actinotignum_timonense", "Adlercreutzia_equolifaciens", "Aeriscardovia_aeriphila",
    "Aerococcus_christensenii", "Aeromonas_caviae", "Agathobaculum_butyriciproducens", "Aggregatibacter_segnis", "Akkermansia_muciniphila",
    "Alistipes_finegoldii", "Alistipes_indistinctus", "Alistipes_inops", "Alistipes_onderdonkii", "Alistipes_putredinis",
    "Alistipes_shahii", "Alistipes_timonensis", "Allisonella_histaminiformans", "Alloprevotella_tannerae", "Alloscardovia_omnicolens",
    "Anaerobiospirillum_thomasii", "Anaerococcus_lactolyticus", "Anaerococcus_vaginalis", "Anaerofilum_sp_An201", "Anaerofustis_stercorihominis",
    "Anaeroglobus_geminatus", "Anaeromassilibacillus_sp_An250", "Anaerosporobacter_mobilis", "Anaerostipes_caccae", "Anaerostipes_hadrus",
    "Anaerostipes_sp_494a", "Anaerotignum_lactatifermentans", "Anaerotruncus_colihominis", "Anaerotruncus_sp_CAG_528", "Asaccharobacter_celatus",
    "Atlantibacter_hermannii", "Atopobium_deltae", "Atopobium_minutum", "Atopobium_parvulum", "Atopobium_rimae", "Atopobium_vaginae", # A - 61 species

    "Bacillus_intestinalis", "Bacteroides_caccae", "Bacteroides_caecimuris", "Bacteroides_cellulosilyticus", "Bacteroides_clarus",
    "Bacteroides_coagulans", "Bacteroides_coprocola", "Bacteroides_coprophilus", "Bacteroides_dorei", "Bacteroides_eggerthii",
    "Bacteroides_faecis", "Bacteroides_faecis_CAG_32", "Bacteroides_finegoldii", "Bacteroides_fluxus", "Bacteroides_fragilis",
    "Bacteroides_galacturonicus", "Bacteroides_intestinalis", "Bacteroides_massiliensis", "Bacteroides_nordii", "Bacteroides_ovatus",
    "Bacteroides_pectinophilus", "Bacteroides_plebeius", "Bacteroides_salyersiae", "Bacteroides_sartorii", "Bacteroides_sp_43_108",
    "Bacteroides_sp_CAG_144", "Bacteroides_sp_CAG_443", "Bacteroides_sp_CAG_462", "Bacteroides_sp_CAG_530", "Bacteroides_sp_CAG_598",
    "Bacteroides_sp_CAG_633", "Bacteroides_sp_D2", "Bacteroides_sp_OM08_11", "Bacteroides_stercorirosoris", "Bacteroides_stercoris",
    "Bacteroides_thetaiotaomicron", "Bacteroides_uniformis", "Bacteroides_vulgatus", "Bacteroides_xylanisolvens", "Bacteroidetes_oral_taxon_274",
    "Barnesiella_intestinihominis", "Barnesiella_sp_An22", "Bifidobacterium_adolescentis", "Bifidobacterium_angulatum", "Bifidobacterium_animalis",
    "Bifidobacterium_bifidum", "Bifidobacterium_breve", "Bifidobacterium_catenulatum", "Bifidobacterium_dentium", "Bifidobacterium_gallinarum",
    "Bifidobacterium_kashiwanohense", "Bifidobacterium_longum", "Bifidobacterium_moukalabense", "Bifidobacterium_pseudocatenulatum", "Bifidobacterium_pseudolongum",
    "Bifidobacterium_pullorum", "Bifidobacterium_ruminantium", "Bifidobacterium_saeculare", "Bifidobacterium_scardovii", "Bilophila_wadsworthia",
    "Blastocystis_sp_subtype_1", "Blautia_coccoides", "Blautia_hansenii", "Blautia_hydrogenotrophica", "Blautia_obeum",
    "Blautia_producta", "Blautia_sp_An249", "Blautia_sp_CAG_257", "Blautia_sp_N6H1_15", "Blautia_wexlerae",
    "Brachyspira_sp_CAG_700", "Bulleidia_extructa", "Butyribacterium_methylotrophicum", "Butyricicoccus_pullicaecorum", "Butyricimonas_synergistica", # B - 78 species
    "Butyricimonas_virosa", "Butyrivibrio_crossotus", "Butyrivibrio_sp_CAG_318",
    
    "Campylobacter_concisus", "Campylobacter_curvus", "Campylobacter_gracilis", "Campylobacter_hominis", "Campylobacter_rectus",
    "Campylobacter_showae", "Campylobacter_upsaliensis", "Campylobacter_ureolyticus", "Candidatus_Gastranaerophilales_bacterium", "Candidatus_Methanomassiliicoccus_intestinalis",
    "Candidatus_Stoquefichus_sp_KLE1796", "Catabacter_hongkongensis", "Catenibacterium_mitsuokai", "Cellulosilyticum_lentocellum", "Christensenella_minuta",
    "Chryseobacterium_arthrosphaerae", "Citrobacter_amalonaticus", "Citrobacter_braakii", "Citrobacter_europaeus", "Citrobacter_farmeri",
    "Citrobacter_freundii", "Citrobacter_koseri", "Citrobacter_pasteurii", "Citrobacter_portucalensis", "Citrobacter_sp_MGH106",
    "Citrobacter_werkmanii", "Citrobacter_youngae", "Cloacibacillus_evryensis", "Cloacibacillus_porcorum", "Clostridiales_bacterium_1_7_47FAA",
    "Clostridiales_bacterium_CHKCI006", "Clostridioides_difficile", "Clostridium_aldenense", "Clostridium_asparagiforme", "Clostridium_bolteae",
    "Clostridium_bolteae_CAG_59", "Clostridium_botulinum", "Clostridium_butyricum", "Clostridium_cadaveris", "Clostridium_celatum",
    "Clostridium_celerecrescens", "Clostridium_citroniae", "Clostridium_clostridioforme", "Clostridium_disporicum", "Clostridium_hiranonis",
    "Clostridium_hylemonae", "Clostridium_innocuum", "Clostridium_lavalense", "Clostridium_leptum", "Clostridium_methylpentosum",
    "Clostridium_neonatale", "Clostridium_paraputrificum", "Clostridium_perfringens", "Clostridium_saccharolyticum", "Clostridium_scindens",
    "Clostridium_sp_7_2_43FAA", "Clostridium_sp_CAG_167", "Clostridium_sp_CAG_242", "Clostridium_sp_CAG_253", "Clostridium_sp_CAG_299",
    "Clostridium_sp_CAG_413", "Clostridium_sp_CAG_58", "Clostridium_sp_CAG_678", "Clostridium_sp_CAG_964", "Clostridium_sp_D5",
    "Clostridium_sp_MSTE9", "Clostridium_sp_chh4_2", "Clostridium_spiroforme", "Clostridium_symbiosum", "Clostridium_ventriculi",
    "Collinsella_aerofaciens", "Collinsella_intestinalis", "Collinsella_massiliensis", "Collinsella_stercoris", "Comamonas_kerstersii",
    "Comamonas_terrigena", "Coprobacillus_cateniformis", "Coprobacter_fastidiosus", "Coprobacter_secundus", "Coprobacter_sp",
    "Coprococcus_catus", "Coprococcus_comes", "Coprococcus_eutactus", "Corynebacterium_amycolatum", "Corynebacterium_argentoratense",
    "Corynebacterium_aurimucosum", "Corynebacterium_variabile", "Cryptococcus_neoformans", "Cutibacterium_acnes", "Cutibacterium_avidum", # C - 90 species

    "Delftia_acidovorans", "Delftia_tsuruhatensis", "Desulfovibrio_fairfieldensis", "Desulfovibrio_legallii", "Desulfovibrio_piger",
    "Desulfovibrionaceae_bacterium", "Dialister_invisus", "Dialister_micraerophilus", "Dialister_pneumosintes", "Dialister_sp_CAG_357",
    "Dialister_succinatiphilus", "Dielma_fastidiosa", "Dorea_formicigenerans", "Dorea_longicatena", "Dorea_sp_CAG_317",
    "Dorea_sp_D27", "Dysgonomonas_sp_37_18", # D - 17 species
    
    "Eggerthella_lenta", "Eikenella_corrodens", "Eisenbergiella_massiliensis", "Eisenbergiella_tayi", "Elizabethkingia_bruuniana",
    "Elizabethkingia_miricola", "Enorma_massiliensis", "Enterobacter_bugandensis", "Enterobacter_cloacae_complex", "Enterobacter_mori",
    "Enterococcus_avium", "Enterococcus_canintestini", "Enterococcus_casseliflavus", "Enterococcus_dispar", "Enterococcus_durans",
    "Enterococcus_faecalis", "Enterococcus_faecium", "Enterococcus_gallinarum", "Enterococcus_gilvus", "Enterococcus_hirae",
    "Enterococcus_italicus", "Enterococcus_malodoratus", "Enterococcus_raffinosus", "Enterococcus_sp_3H8_DIV0648", "Enterorhabdus_caecimuris",
    "Erysipelatoclostridium_ramosum", "Erysipelothrix_larvae", "Escherichia_albertii", "Escherichia_coli", "Escherichia_marmotae",
    "Eubacterium_brachy", "Eubacterium_callanderi", "Eubacterium_coprostanoligenes", "Eubacterium_dolichum_CAG_375", "Eubacterium_eligens",
    "Eubacterium_hallii", "Eubacterium_limosum", "Eubacterium_ramulus", "Eubacterium_rectale", "Eubacterium_siraeum",
    "Eubacterium_sp_An11", "Eubacterium_sp_CAG_180", "Eubacterium_sp_CAG_251", "Eubacterium_sp_CAG_274", "Eubacterium_sp_CAG_38",
    "Eubacterium_sp_OM08_24", "Eubacterium_sulci", "Eubacterium_ventriosum", # E - 48 species
    
    "Faecalibacterium_prausnitzii", "Faecalicatena_orotica", "Faecalicoccus_pleomorphus", "Faecalitalea_cylindroides", "Finegoldia_magna",
    "Firmicutes_bacterium_CAG_110", "Firmicutes_bacterium_CAG_145", "Firmicutes_bacterium_CAG_170", "Firmicutes_bacterium_CAG_424", "Firmicutes_bacterium_CAG_646",
    "Firmicutes_bacterium_CAG_791", "Firmicutes_bacterium_CAG_83", "Firmicutes_bacterium_CAG_94", "Firmicutes_bacterium_CAG_95", "Flavonifractor_plautii",
    "Flavonifractor_sp_An10", "Flavonifractor_sp_An100", "Fusicatenibacter_saccharivorans", "Fusobacterium_gonidiaformans", "Fusobacterium_hwasookii",
    "Fusobacterium_mortiferum", "Fusobacterium_naviforme", "Fusobacterium_necrophorum", "Fusobacterium_nucleatum", "Fusobacterium_periodonticum",
    "Fusobacterium_sp_CAG_439", "Fusobacterium_sp_oral_taxon_370", "Fusobacterium_ulcerans", "Fusobacterium_varium", # F - 29 species
    
    "Gardnerella_vaginalis", "Gemella_asaccharolytica", "Gemella_haemolysans", "Gemella_morbillorum", "Gemella_sanguinis",
    "Gemmiger_formicilis", "Gemmiger_sp_An50", "Gemmiger_sp_An87", "Gordonia_jacobaea", "Gordonia_sputi",
    "Gordonibacter_pamelaeae", "Granulicatella_elegans", # G = 12 species
    
    "Haemophilus_haemolyticus", "Haemophilus_parahaemolyticus", "Haemophilus_parainfluenzae", "Haemophilus_paraphrohaemolyticus", "Haemophilus_pittmaniae",
    "Haemophilus_sp_HMSC71H05", "Haemophilus_sputorum", "Hafnia_alvei", "Hafnia_paralvei", "Harryflintia_acetispora",
    "Helicobacter_bilis", "Helicobacter_canis", "Holdemanella_biformis", "Holdemania_filiformis", "Hungatella_effluvii",
    "Hungatella_hathewayi", # H - 16 species
    
    "Intestinibacter_bartlettii", "Intestinimonas_butyriciproducens", # I - 2 species
    
    "Klebsiella_aerogenes", "Klebsiella_grimontii", "Klebsiella_michiganensis", "Klebsiella_oxytoca", "Klebsiella_pneumoniae",
    "Klebsiella_quasipneumoniae", "Klebsiella_quasivariicola", "Klebsiella_variicola", "Kluyvera_ascorbata", "Kluyvera_cryocrescens",
    "Kluyvera_georgiana", "Kluyvera_intestini", "Kosakonia_sacchari", "Kosakonia_sp_S29", # K - 14 species
    
    "Lachnoclostridium_sp_An118", "Lachnoclostridium_sp_An131", "Lachnoclostridium_sp_An138", "Lachnoclostridium_sp_An14", "Lachnoclostridium_sp_An169",
    "Lachnoclostridium_sp_An181", "Lachnoclostridium_sp_An76", "Lachnospira_pectinoschiza", "Lachnospiraceae_bacterium_2_1_46FAA", "Lachnospiraceae_bacterium_oral_taxon_096",
    "Lactobacillus_acidophilus", "Lactobacillus_algidus", "Lactobacillus_amylovorus", "Lactobacillus_animalis", "Lactobacillus_buchneri",
    "Lactobacillus_crispatus", "Lactobacillus_curvatus", "Lactobacillus_delbrueckii", "Lactobacillus_fermentum", "Lactobacillus_gallinarum",
    "Lactobacillus_gasseri", "Lactobacillus_helveticus", "Lactobacillus_iners", "Lactobacillus_johnsonii", "Lactobacillus_kalixensis",
    "Lactobacillus_mucosae", "Lactobacillus_oris", "Lactobacillus_paragasseri", "Lactobacillus_plantarum", "Lactobacillus_reuteri",
    "Lactobacillus_rhamnosus", "Lactobacillus_rogosae", "Lactobacillus_ruminis", "Lactobacillus_sakei", "Lactobacillus_salivarius",
    "Lactobacillus_vaginalis", "Lactococcus_garvieae", "Lactococcus_lactis", "Lactococcus_petauri", "Lactococcus_raffinolactis",
    "Lactonifactor_longoviformis", "Lawsonella_clevelandensis", "Lawsonibacter_asaccharolyticus", "Leclercia_adecarboxylata", "Leuconostoc_garlicum",
    "Leuconostoc_gelidum", "Leuconostoc_lactis", "Leuconostoc_mesenteroides", "Leuconostoc_pseudomesenteroides", "Listeria_monocytogenes", # L - 50 species

    "Mageeibacillus_indolicus", "Massiliomicrobiota_timonensis", "Megamonas_funiformis", "Megamonas_funiformis_CAG_377", "Megamonas_hypermegale",
    "Megasphaera_elsdenii", "Megasphaera_hexanoica", "Megasphaera_micronuciformis", "Megasphaera_sp_BV3C16_1", "Megasphaera_sp_DISK_18",
    "Megasphaera_sp_MJR8396C", "Methanobrevibacter_smithii", "Methanosphaera_stadtmanae", "Mitsuokella_jalaludinii", "Mitsuokella_multacida",
    "Mogibacterium_diversum", "Mogibacterium_pumilum", "Monoglobus_pectinilyticus", "Morganella_morganii", "Muribaculum_intestinale",
    "Murimonas_intestini", "Mycoplasma_hominis", "Negativicoccus_succinicivorans", "Neisseria_flavescens", "Neisseria_subflava", # M - 25 species
    
    "Odoribacter_laneus", "Odoribacter_splanchnicus", "Olsenella_scatoligenes", "Oscillibacter_sp_57_20", "Oscillibacter_sp_CAG_241",
    "Oxalobacter_formigenes", # O - 6 species
    
    "Paeniclostridium_sordellii", "Pantoea_agglomerans", "Pantoea_dispersa", "Parabacteroides_distasonis", "Parabacteroides_goldsteinii",
    "Parabacteroides_gordonii", "Parabacteroides_johnsonii", "Parabacteroides_merdae", "Parabacteroides_sp_CAG_409", "Paraprevotella_clara",
    "Paraprevotella_xylaniphila", "Parasutterella_excrementihominis", "Parvimonas_micra", "Parvimonas_sp_KA00067", "Parvimonas_sp_oral_taxon_110",
    "Pediococcus_acidilactici", "Peptococcus_niger", "Peptoniphilus_coxii", "Peptoniphilus_duerdenii", "Peptoniphilus_harei",
    "Peptoniphilus_lacrimalis", "Peptoniphilus_sp_BV3C26", "Peptoniphilus_sp_HMSC062D09", "Peptostreptococcus_anaerobius", "Peptostreptococcus_sp_MV1",
    "Peptostreptococcus_stomatis", "Phascolarctobacterium_faecium", "Phascolarctobacterium_sp_CAG_266", "Phascolarctobacterium_succinatutens", "Phytobacter_ursingii",
    "Plesiomonas_shigelloides", "Pluralibacter_gergoviae", "Porphyromonas_asaccharolytica", "Porphyromonas_endodontalis", "Porphyromonas_gingivalis",
    "Porphyromonas_sp_HMSC065F10", "Porphyromonas_uenonis", "Prevotella_amnii", "Prevotella_bergensis", "Prevotella_bivia",
    "Prevotella_buccae", "Prevotella_buccalis", "Prevotella_colorans", "Prevotella_copri", "Prevotella_corporis",
    "Prevotella_denticola", "Prevotella_disiens", "Prevotella_histicola", "Prevotella_intermedia", "Prevotella_melaninogenica",
    "Prevotella_nigrescens", "Prevotella_salivae", "Prevotella_sp_885", "Prevotella_sp_AM42_24", "Prevotella_sp_CAG_1031",
    "Prevotella_sp_CAG_1092", "Prevotella_sp_CAG_1185", "Prevotella_sp_CAG_1320", "Prevotella_sp_CAG_279", "Prevotella_sp_CAG_520",
    "Prevotella_sp_CAG_5226", "Prevotella_sp_CAG_873", "Prevotella_sp_CAG_891", "Prevotella_sp_oral_taxon_306", "Prevotella_stercorea",
    "Prevotella_timonensis", "Propionibacterium_freudenreichii", "Proteobacteria_bacterium_CAG_139", "Proteus_hauseri", "Proteus_mirabilis",
    "Providencia_rustigianii", "Pseudescherichia_vulneris", "Pseudoflavonifractor_capillosus", "Pseudoflavonifractor_sp_An184", "Pseudomonas_aeruginosa_group",
    "Pyramidobacter_piscolens", "Pyramidobacter_sp_C12_8", # P - 77 species
    
    "Raoultella_ornithinolytica", "Raoultella_planticola", "Rikenella_microfusus", "Robinsoniella_sp_RHS", "Romboutsia_ilealis",
    "Roseburia_faecis", "Roseburia_hominis", "Roseburia_intestinalis", "Roseburia_inulinivorans", "Roseburia_sp_CAG_182",
    "Roseburia_sp_CAG_303", "Roseburia_sp_CAG_309", "Roseburia_sp_CAG_471", "Rothia_mucilaginosa", "Ruminococcaceae_bacterium_D16",
    "Ruminococcaceae_bacterium_D5", "Ruminococcus_bicirculans", "Ruminococcus_bromii", "Ruminococcus_callidus", "Ruminococcus_champanellensis",
    "Ruminococcus_gnavus", "Ruminococcus_lactaris", "Ruminococcus_obeum_CAG_39", "Ruminococcus_sp_CAG_330", "Ruminococcus_sp_CAG_403",
    "Ruminococcus_sp_CAG_488", "Ruminococcus_sp_CAG_624", "Ruminococcus_torques", "Ruthenibacterium_lactatiformans", # R - 29 species
    
    "Saccharomyces_cerevisiae", "Salmonella_enterica", "Sanguibacteroides_justesenii", "Scardovia_wiggsiae", "Sellimonas_intestinalis",
    "Serratia_liquefaciens", "Serratia_marcescens", "Slackia_isoflavoniconvertens", "Sneathia_amnii", "Solobacterium_moorei",
    "Staphylococcus_argenteus", "Staphylococcus_aureus", "Staphylococcus_epidermidis", "Staphylococcus_hominis", "Staphylococcus_lugdunensis",
    "Stenotrophomonas_maltophilia", "Streptococcus_agalactiae", "Streptococcus_anginosus_group", "Streptococcus_australis", "Streptococcus_gallolyticus",
    "Streptococcus_gordonii", "Streptococcus_infantarius", "Streptococcus_infantis", "Streptococcus_lutetiensis", "Streptococcus_milleri",
    "Streptococcus_mitis", "Streptococcus_mutans", "Streptococcus_oralis", "Streptococcus_parasanguinis", "Streptococcus_parauberis",
    "Streptococcus_pasteurianus", "Streptococcus_peroris", "Streptococcus_pneumoniae", "Streptococcus_pseudopneumoniae", "Streptococcus_salivarius",
    "Streptococcus_salivarius_CAG_79", "Streptococcus_sanguinis", "Streptococcus_sobrinus", "Streptococcus_sp_A12", "Streptococcus_sp_F0442",
    "Streptococcus_sp_HPH0090", "Streptococcus_sp_M334", "Streptococcus_sp_SK643", "Streptococcus_thermophilus", "Streptococcus_vestibularis",
    "Succinatimonas_hippei", "Sutterella_parvirubra", # S - 47 species
    
    "Tannerella_forsythia", "Terrisporobacter_othiniensis", "Tissierellia_bacterium_KA00581", "Treponema_succinifaciens", "Trueperella_bernardiae",
    "Turicibacter_sanguinis", "Turicimonas_muris", "Tyzzerella_nexilis", # T - 8 species
    
    "Ureaplasma_parvum", # U - 1 species
    
    "Varibaculum_cambriense", "Veillonella_atypica", "Veillonella_dispar", "Veillonella_infantium", "Veillonella_parvula",
    "Veillonella_rodentium", "Veillonella_rogosae", "Veillonella_seminalis", "Veillonella_sp_CAG_933", "Veillonella_sp_T11011_6",
    "Veillonella_tobetsuensis", "Veillonellaceae_bacterium_DNF00626", "Victivallis_vadensis", # V - 13 species
    
    "Weissella_cibaria", "Weissella_confusa", # W - 2 species
]

function fix_metadata_colnames!(longdata_df::DataFrame)

    ## This function was created to fix small inconsistencies and typos on the longdata table.
    ## As original data cleanup and wrangling progresses, it will be altered and ultimately deprecated.
    ## Comments are provided to help review transformations

    # 1. Rename "[Collinsella]_massiliensis" to "Collinsella_massiliensis"
    longdata_df[longdata_df.variable .== "[Collinsella]_massiliensis", :variable] .= "Collinsella_massiliensis"

    return(longdata_df)
    
end

function check_longdata_metaduplicates!(longdata_df::DataFrame; remove_duplicates=true)

    n_unique_rows = size(unique(longdata_df[:, 1:3]),1)
    
    if (n_unique_rows != size(longdata_df, 1))
    
        n_nonunique = size(longdata_df, 1) - n_unique_rows

        if remove_duplicates

            @warn "Long Dataframe contains non-unique rows! Argument `remove_duplicates` set to `true`. Removing duplicated data."
            @warn "After removal, $(n_unique_rows) will remain. $(n_nonunique) rows were rmoved from the original $(size(longdata_df, 1))"
            rawDf = unique!(longdata_df)
            return(longdata_df)

        else

            @warn "Long Dataframe contains non-unique rows! Argument `remove_duplicates` set to `false`. Returning source data."
            return(longdata_df)

        end #end if remove_duplicates
    
    end # end if n_unique_rows
     
end # end function

function build_future_df(base_df, to_predict::Symbol)

    subjects = unique(base_df.subject)

    result_lines = Vector{DataFrame}()

    for this_subject in subjects

        subject_df = base_df[ base_df.subject .== this_subject ,:]

        if (size(subject_df, 1) == 1)

            continue;

        else

            for origin_idx in 1:(size(subject_df, 1) - 1)

                if (ismissing(subject_df[origin_idx, :Absiella_dolichum]))

                    continue;

                else

                    for target_idx in (origin_idx + 1):size(subject_df, 1)

                        if (ismissing(subject_df[target_idx, to_predict]))

                            continue;

                        else
                            origin_df = DataFrame()
                            push!(origin_df, copy(subject_df[origin_idx, :]))
                            target_df = subject_df[target_idx, :]

                            insertcols!(origin_df, 1, :target => parse(Float64, target_df[to_predict]))
                            insertcols!(origin_df, 1, :futureAgeMonths => parse(Float64, target_df[:ageMonths]))
                            insertcols!(origin_df, 1, :ageMonthsDelta => parse(Float64, target_df[:ageMonths]) - parse(Float64, origin_df[1, :ageMonths]))
                            insertcols!(origin_df, 1, :futureTimepoint => subject_df[target_idx, :timepoint])
                            insertcols!(origin_df, 1, :timepointDelta => subject_df[target_idx, :timepoint] - subject_df[origin_idx, :timepoint])

                            push!(result_lines, origin_df)

                        end # end if ismissing(to_predict from irigin_idx)

                    end # end for target_idx

                end # end if ismissing(stool from irigin_idx)

            end # end for origin_idx

        end # end if size == 1

    end # end for subject

    return(reduce(vcat, result_lines))

end # end function1