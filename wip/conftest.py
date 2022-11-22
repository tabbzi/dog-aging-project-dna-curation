import pytest
from io import StringIO
from contextlib import contextmanager


@pytest.fixture()
def variants_of_interest():
    return StringIO(
        'tab,trait,phenotype,locus,variant,module,gene,title,desc,CHR,POS,normalAllele,variantAllele,imputation,reference\n'
        'CoatColor,black,liver,B,bc,mod_liver,TYRP1,Liver - variant p.(C41S),Alters eumelanin pigment to a brown hue,11,33317810,T,A,impute-v2,CanFam3.1\n'
        'CoatColor,black,liver,B,bd,mod_liver,TYRP1,Liver - variant p.(P345del),Alters eumelanin pigment to a brown hue,11,33326726,ACCT,A,impute-v2,CanFam3.1\n'
        'CoatColor,black,liver,B,bs,mod_liver,TYRP1,Liver - variant p.(Gln331*),Alters eumelanin pigment to a brown hue,11,33326685,C,T,impute-v2,CanFam3.1\n'
        'CoatColor,black,cocoa,Co,co,mod_cocoa,HPS3,Cocoa - variant p.(T807*),Alters eumelanin pigment to a deep brown hue,23,43969695,G,A,impute-v2,CanFam3.1\n'
        'CoatColor,black,dilution,D,d1,mod_dilute,MLPH,Dilution - splice variant,Dilutes eumelanin pigment to a lighter hue,25,48121642,G,A,impute-v2,CanFam3.1\n'
        'CoatColor,black,dilution,D,d2,mod_dilute,MLPH,Dilution - variant p.(Q235H),Dilutes eumelanin pigment to a lighter hue,25,48150787,G,C,impute-v2,CanFam3.1\n'
        'CoatColor,black,dilution,D,d3,mod_dilute,MLPH,Dilution - variant p.(H223Pfs*41),Dilutes eumelanin pigment to a lighter hue,25,48150749,C,CC,unobserved,CanFam3.1\n'
        'CoatColor,red,red_intensity,I,i1,mod_red_intensity,lincRNA,Red intensity - marker 1,Influences deepness of pheomelanin pigment,2,74746906,T,A,impute-v2,CanFam3.1\n'
        'CoatColor,red,red_intensity,I,i2,mod_red_intensity,intergenic,Red intensity - marker 2,Influences deepness of pheomelanin pigment,15,29840789,T,C,impute-v2,CanFam3.1\n'
        'CoatColor,red,red_intensity,I,i3,mod_red_intensity,SLC264A,Red intensity - marker 3,Influences deepness of pheomelanin pigment,18,12910382,T,C,impute-v2,CanFam3.1\n'
        'CoatColor,red,red_intensity,I,i4,mod_red_intensity,intergenic,Red intensity - marker 4,Influences deepness of pheomelanin pigment,20,55850145,T,C,impute-v2,CanFam3.1\n'
        'CoatColor,red,red_intensity,I,i5,mod_red_intensity,TYR,Red intensity - marker 5,Influences deepness of pheomelanin pigment,21,10864834,G,A,impute-v2,CanFam3.1\n'
        'CoatPattern,agouti,sable,A,Ay1,mod_sable,ASIP,Sable - variant p.(A82S),Determines if a dog has sable pattern,24,23393510,G,T,impute-v2,CanFam3.1\n'
        'CoatPattern,agouti,sable,A,Ay2,mod_sable,ASIP,Sable - variant p.(R83H),Determines if a dog has sable pattern,24,23393514,G,A,impute-v2,CanFam3.1\n'
        'CoatPattern,agouti,tan_points,A,at,mod_tan_points,ASIP,Tan points - marker,Determines if a dog has tan point pattern,24,23365205,C,T,impute-v2,CanFam3.1\n'
        'CoatPattern,agouti,recessive_black,A,a-,mod_rec_black,ASIP,Recessive black - variant p.(R96C),Eumelanin solidly overrides pheomelanin across entire coat,24,23393552,C,T,impute-v2,CanFam3.1\n'
        'CoatPattern,agouti,saddle,A,atm1,mod_saddle,RALY,Saddle - marker 1,Allows for saddle pattern to supercede tan points,24,23252763,CAGAGTTTCCCCAGGT,C,impute-v2,CanFam3.1\n'
        'CoatPattern,agouti,saddle,A,atm2,mod_saddle,RALY,Saddle - marker 2,Allows for saddle pattern to supercede tan point pattern,24,23252753,GTCCCCAGGTCAGAGTT,G,impute-v2,CanFam3.1\n'
        'CoatPattern,mask,melanistic_mask,E,EM,mod_mask,MC1R,Facial mask - variant p.(M264V),Determines if a dog has a dark facial mask,5,63694460,T,C,impute-v2,CanFam3.1\n'
        'CoatPattern,grizzle,sighthound_domino,E,Eg,mod_grizzle,MC1R,Sighthound grizzle - variant p.(G78V),Determines if dog has shaded sable with widows peak pattern,5,63695017,C,A,impute-v2,CanFam3.1\n'
        'CoatPattern,northern_domino,northern_domino,E,EA,mod_domino,MC1R,Northern domino - variant p.(R301C),Determines if a dog has domino facial and body pattern common in Northern breeds,5,63694349,G,A,impute-v2,CanFam3.1\n'
        'CoatPattern,red,recessive_red,E,e1,mod_rec_red,MC1R,Recessive red - variant p.(R306*),Pheomelanin solidly overrides eumelanin across entire coat,5,63694334,G,A,impute-v2,CanFam3.1\n'
        'CoatPattern,red,recessive_red,E,e2,mod_rec_red,MC1R,Recessive red - regulatory variant,Pheomelanin solidly overrides eumelanin across entire coat,5,63695679,C,G,impute-v2,CanFam3.1\n'
        'CoatPattern,brindle,dominant_black,K,KB,mod_dom_black,CBD103,Dominant black - variant p.(G23del),Eumelanin solidly overrides pheomelanin across entire coat,16,58965448,TCCC,T,impute-v2,CanFam3.1\n'
        'CoatPattern,brindle,brindle,K,Kb1,mod_brindle,intergenic,Brindle - marker 1,Associated with brindle coat pattern,16,58999689,A,AGG,impute-v2,CanFam3.1\n'
        'CoatPattern,brindle,brindle,K,Kb2,mod_brindle,intergenic,Brindle - marker 2,Associated with brindle coat pattern,16,58999678,GCTTCCCTAAAA,G,impute-v2,CanFam3.1\n'
        'CoatPattern,roan,ticking,T,ta,mod_ticking,USH2A,Ticking - marker,,38,11165134,G,A,impute-v2,CanFam3.1\n'
        'CoatPattern,merle,harlequin,H,h,mod_harlequin,PSMB7,Harlequin - variant p.(V49I),Increases the amount of white on a merle patterned dog,9,58530295,T,G,impute-v2,CanFam3.1\n'
        'CoatType,coat_texture,curly_coat,Cu,c1,mod_curly_coat,KRT71,Curly coat - variant p.(R151W),Makes coat curly or wavy,27,2539211,C,T,impute-v2,CanFam3.1\n'
        'CoatType,coat_texture,curly_coat,Cu,c2,mod_curly_coat,KRT71,Curly coat - variant p.(S422Rfs),Makes coat curly or wavy,27,2543231,CTG,C,impute-v2,CanFam3.1\n'
        'CoatType,coat_length,long_coat,L,Lh,mod_coat_length,FGF5,Long coat - variant p.(C95F),Increases coat length,32,4509367,G,T,impute-v2,CanFam3.1\n'
        'CoatType,furnishings,furnishings,IC,ins,mod_furnishings,RSPO2,Furnishings - marker,Causes long hair on eyebrows and jowls,13,8491477,A,C,impute-v2,CanFam3.1\n'
        'CoatType,shedding,shedding_propensity,sp,sp,Shedding,MC5R,Shedding propensity - variant p.(A237T),Decreases shedding,1,24430748,T,C,impute-v2,CanFam3.1\n'
        'CoatType,coat_layer,single_coat,sc,1,mod_coat_layer,ADRB1-AU1,Single-layer coat - marker 1,Loss of double-layered coat,28,24860187,C,T,impute-v2,CanFam3.1\n'
        'CoatType,coat_layer,single_coat,sc,2,mod_coat_layer,ADRB1-AU1,Single-layer coat - marker 2,Loss of double-layered coat,28,24870184,G,A,impute-v2,CanFam3.1\n'
        'SpecialFeatures,altitude,high_altitude_adaptation,alt1,a1,AltitudeAdaptation,EPAS1,High altitude hypoxia tolerance - marker 1,"Observed in breeds from high elevations, like Tibetan mastiffs",10,48626862,G,A,impute-v2,CanFam3.1\n'
        'SpecialFeatures,altitude,high_altitude_adaptation,alt2,a2,AltitudeAdaptation,EPAS1,High altitude hypoxia tolerance - marker 2,"Observed in breeds from high elevations, like Tibetan mastiffs",10,48630137,G,T,impute-v2,CanFam3.1\n'
        'SpecialFeatures,altitude,high_altitude_adaptation,alt3,a3,AltitudeAdaptation,EPAS1,High altitude hypoxia tolerance - marker 3,"Observed in breeds from high elevations, like Tibetan mastiffs",10,48630153,G,A,impute-v2,CanFam3.1\n'
        'SpecialFeatures,altitude,high_altitude_adaptation,alt4,a4,AltitudeAdaptation,EPAS1,High altitude hypoxia tolerance - marker 4,"Observed in breeds from high elevations, like Tibetan mastiffs",10,48633379,C,T,impute-v2,CanFam3.1\n'
        'SpecialFeatures,snout_length,brachycephaly,bmp3,br1,Brachycephaly,BMP3,Brachycephaly - variant p.(F452L),"Observed in breeds with short snouts, like Pugs",32,5231894,C,A,impute-v2_disabled,CanFam3.1\n'
        'SpecialFeatures,eye_color,blue_eyes,alx4,blue,mod_eye_color,ALX4,Blue eyes - marker,Associated with blue-colored eyes,18,44924848,C,T,impute-v2,CanFam3.1\n'
        'SpecialFeatures,food_motivation,food_motivation,POMC,del,FoodMotivation,POMC,Food motivation - marker,Observed in Labradors predisposed to obesity,17,19431807,CGCGGCGGGGCCCT,C,impute-v2_disabled,CanFam3.1\n'
        'SpecialFeatures,leg_length,short_legs,fgf4c18gwas,snp,mod_short_legs,FGF4 retrogene on chromosome 18,Shortened legs - marker,"Observed in breeds with short legs, such as Corgis",18,20423056,A,G,impute-v2,CanFam3.1\n'
        'SpecialFeatures,leg_length,long_legs,esr1,e1,mod_long_legs,ESR1,Long legs - marker 1,Associated with longer than average legs,1,42302476,C,T,impute-v2,CanFam3.1\n'
        'SpecialFeatures,leg_length,long_legs,esr1,e2,mod_long_legs,ESR1,Long legs - marker 2,Associated with longer than average legs,1,42306135,A,G,impute-v2,CanFam3.1\n'
        'SpecialFeatures,leg_length,long_legs,esr1,e3,mod_long_legs,ESR1,Long legs - marker 3,Associated with longer than average legs,1,42300805,T,G,impute-v2,CanFam3.1\n'
        'SpecialFeatures,tail_length,natural_bob_tail,BT,bt,BobTail,T,Natural bob tail - variant p.(I63M),Dogs with variant born with shorter tail,1,54192143,G,C,impute-v2,CanFam3.1\n'
        'wrong imputation,tail_length,natural_bob_tail,BT,bt,BobTail,T,Natural bob tail - variant p.(I63M),Dogs with variant born with shorter tail,1,54192143,G,C,impute-v1,CanFam3.1\n'
        'wrong ref,tail_length,natural_bob_tail,BT,bt,BobTail,T,Natural bob tail - variant p.(I63M),Dogs with variant born with shorter tail,1,54192143,G,C,impute-v2,CanFam4.1\n'
    )

@pytest.fixture()
def sample_vcf():
    return StringIO(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t1\t2\t3\t4\t5\n'
        '1\t24430748\t1:24430748:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t0/0\t0/1\t1/1\t0/0\t0/1\n'
        '1\t42300796\t1:42300796:G:GT\tG\tGT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42300804\t1:42300804:TG:T\tTG\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42300805\t1:42300805:G:T\tG\tT\t.\tLOWCONF\t.\tGT\t0/0\t1/0\t1/0\t0/0\t0/0\n'
        '1\t42300805\t1:42300805:G:*\tG\t*\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0\t0/0\n'
        '1\t42302476\t1:42302476:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42302488\t1:42302488:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42306118\t1:42306118:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42306122\t1:42306122:G:T\tG\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42306135\t1:42306135:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42306155\t1:42306155:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t42306170\t1:42306170:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t54192094\t1:54192094:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t54192102\t1:54192102:TC:T\tTC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '1\t54192143\t1:54192143:G:C\tG\tC\t.\tPASS\t.\tGT\t./.\t0/0\t0/0\t0/0\t0/0\n'
        '2\t74746888\t2:74746888:G:C\tG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '2\t74746891\t2:74746891:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '2\t74746906\t2:74746906:T:A\tT\tA\t.\tLOWCONF\t.\tGT\t0/1\t1/1\t0/1\t1/1\t0/0\n'
        '2\t74746943\t2:74746943:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694334\t5:63694334:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/1\t1/1\t0/0\t0/1\t0/0\n'
        '5\t63694334\t5:63694334:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/1\t1/1\t0/0\t0/1\t0/0\n'
        '5\t63694349\t5:63694349:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694349\t5:63694349:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694382\t5:63694382:TGCAGATGATGAGGGTGA:T\tTGCAGATGATGAGGGTGA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694382\t5:63694382:TGCAGATGATGAGGGTGA:T\tTGCAGATGATGAGGGTGA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694418\t5:63694418:G:T\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694432\t5:63694432:CAG:C\tCAG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694442\t5:63694442:G:C\tG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63694460\t5:63694460:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/1\t1/1\t0/0\t1/1\t0/1\n'
        '5\t63694982\t5:63694982:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t0/0\t0/1\t0/0\n'
        '5\t63695000\t5:63695000:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63695001\t5:63695001:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63695017\t5:63695017:C:A\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63695653\t5:63695653:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63695679\t5:63695679:C:G\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63695696\t5:63695696:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '5\t63695697\t5:63695697:G:T\tG\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '9\t58530295\t9:58530295:T:G\tT\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '9\t58530314\t9:58530314:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '9\t58530341\t9:58530341:G:T\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '9\t58530342\t9:58530342:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '9\t58530345\t9:58530345:G:T\tG\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\n'
        '10\t48626846\t10:48626846:A:G\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626862\t10:48626862:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626879\t10:48626879:T:G\tT\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626882\t10:48626882:A:G\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626885\t10:48626885:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626888\t10:48626888:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626891\t10:48626891:G:C\tG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626897\t10:48626897:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626906\t10:48626906:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48626912\t10:48626912:C:A\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630092\t10:48630092:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630137\t10:48630137:T:G\tT\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630137\t10:48630137:T:G\tT\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630153\t10:48630153:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630153\t10:48630153:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630158\t10:48630158:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48630158\t10:48630158:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48633358\t10:48633358:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48633379\t10:48633379:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '10\t48633417\t10:48633417:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '11\t33317810\t11:33317810:T:A\tT\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '11\t33326685\t11:33326685:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\t0/0\n'
        '11\t33326685\t11:33326685:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\t0/0\n'
        '11\t33326726\t11:33326726:ACCT:A\tACCT\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '11\t33326726\t11:33326726:ACCT:A\tACCT\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '11\t33326760\t11:33326760:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '13\t8491466\t13:8491466:T:C\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t0/1\t0/0\t0/1\n'
        '13\t8491477\t13:8491477:A:C\tA\tC\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '13\t8491506\t13:8491506:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t0/1\t0/1\t0/1\t0/0\t0/1\n'
        '15\t29840760\t15:29840760:CTTCT:C\tCTTCT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '15\t29840789\t15:29840789:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t1/1\t0/1\t0/1\t0/1\t1/1\n'
        '16\t58965417\t16:58965417:T:G\tT\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '16\t58965448\t16:58965448:TCCC:T\tTCCC\tT\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t0/1\t0/1\t0/0\n'
        '16\t58965479\t16:58965479:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '16\t58965491\t16:58965491:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '16\t58999678\t16:58999678:GCTTCCCTAAAA:G\tGCTTCCCTAAAA\tG\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '16\t58999678\t16:58999678:GCTTCCCTAAAA:G\tGCTTCCCTAAAA\tG\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '16\t58999689\t16:58999689:A:AGG\tA\tAGG\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '16\t58999689\t16:58999689:A:AGG\tA\tAGG\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '17\t19431777\t17:19431777:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '17\t19431807\t17:19431807:CGCGGCGGGGCCCT:C\tCGCGGCGGGGCCCT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '17\t19431813\t17:19431813:G:C\tG\tC\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1\t1/1\n'
        '17\t19431813\t17:19431813:G:*\tG\t*\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '18\t12910382\t18:12910382:T:C\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\n'
        '18\t20423055\t18:20423055:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\n'
        '18\t20423056\t18:20423056:A:G\tA\tG\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t1/1\t0/0\n'
        '18\t20423080\t18:20423080:A:G\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '18\t44924848\t18:44924848:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0\t0/0\n'
        '18\t44924889\t18:44924889:T:G\tT\tG\t.\tLOWCONF\t.\tGT\t0/1\t1/1\t0/0\t0/0\t0/1\n'
        '18\t44924893\t18:44924893:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/1\t1/1\t0/0\t0/0\t0/1\n'
        '20\t55850141\t20:55850141:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '20\t55850145\t20:55850145:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '21\t10864803\t21:10864803:A:T\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '21\t10864834\t21:10864834:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/1\t0/1\t1/1\t0/0\n'
        '23\t43969675\t23:43969675:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '23\t43969695\t23:43969695:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '24\t23252753\t24:23252753:G:GTCCCCAGGTCAGAGTT\tG\tGTCCCCAGGTCAGAGTT\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t0/0\t0/1\t1/1\n'
        '24\t23252753\t24:23252753:G:GTCCCCAGGTCAGAGTT\tG\tGTCCCCAGGTCAGAGTT\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t0/0\t0/1\t1/1\n'
        '24\t23252763\t24:23252763:C:CAGAGTTTCCCCAGGT\tC\tCAGAGTTTCCCCAGGT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '24\t23252763\t24:23252763:C:CAGAGTTTCCCCAGGT\tC\tCAGAGTTTCCCCAGGT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '24\t23252810\t24:23252810:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '24\t23365205\t24:23365205:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393510\t24:23393510:T:G\tT\tG\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393510\t24:23393510:T:G\tT\tG\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393510\t24:23393510:T:G\tT\tG\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393514\t24:23393514:A:G\tA\tG\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393514\t24:23393514:A:G\tA\tG\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393514\t24:23393514:A:G\tA\tG\t.\tLOWCONF\t.\tGT\t0/0\t1/1\t1/1\t0/1\t1/1\n'
        '24\t23393552\t24:23393552:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/0\n'
        '24\t23393552\t24:23393552:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/0\n'
        '24\t23393552\t24:23393552:C:T\tC\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/0\n'
        '24\t23393575\t24:23393575:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48121642\t25:48121642:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48121657\t25:48121657:C:A\tC\tA\t.\tPASS\t.\tGT\t1/1\t0/0\t1/1\t1/1\t0/1\n'
        '25\t48150710\t25:48150710:G:T\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150741\t25:48150741:A:AC\tA\tAC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150741\t25:48150741:A:AC\tA\tAC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150749\t25:48150749:C:A\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150749\t25:48150749:C:A\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150751\t25:48150751:C:T\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150751\t25:48150751:C:T\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150787\t25:48150787:G:C\tG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150787\t25:48150787:G:C\tG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150807\t25:48150807:GC:G\tGC\tG\t.\tLOWCONF\t.\tGT\t0/1\t0/1\t0/0\t0/1\t0/1\n'
        '25\t48150822\t25:48150822:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '25\t48150823\t25:48150823:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/1\t0/1\t0/1\t0/1\n'
        '25\t48150836\t25:48150836:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '27\t2539211\t27:2539211:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/0\n'
        '27\t2543198\t27:2543198:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '27\t2543230\t27:2543230:C:A\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '27\t2543231\t27:2543231:CTG:C\tCTG\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '27\t2543234\t27:2543234:AAGC:A\tAAGC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '27\t2543257\t27:2543257:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '28\t24860187\t28:24860187:T:C\tT\tC\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t1/1\t0/0\t0/1\n'
        '28\t24860210\t28:24860210:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '28\t24860225\t28:24860225:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t0/0\t0/1\n'
        '28\t24870184\t28:24870184:A:G\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1\t0/0\t0/1\n'
        '32\t4509367\t32:4509367:G:T\tG\tT\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t1/1\t0/0\t0/0\n'
        '32\t4509377\t32:4509377:CATCGGT:C\tCATCGGT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t4509403\t32:4509403:C:G\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t5231858\t32:5231858:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t5231861\t32:5231861:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t5231876\t32:5231876:T:C\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t5231882\t32:5231882:T:C\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t5231888\t32:5231888:C:T\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        '32\t5231894\t32:5231894:C:A\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\t0/0\n'
        '38\t11165134\t38:11165134:G:A\tG\tA\t.\tLOWCONF\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/1\n'
        '38\t11165150\t38:11165150:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
        'chr38\t11165151\t38:11165150:G:A\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n'
    )

@pytest.fixture()
def old_predict_output():
    return StringIO(
        'sample,tab,trait,result,string,value,color,image\n'
        '1,CoatColor,black,black,Black,,#17181D,\n'
        '1,CoatColor,red,tan,Tan,,#A86B39,\n'
        '1,CoatPattern,agouti,sable,Sable,,,CoatPattern_locusA_sable.svg\n'
        '1,CoatPattern,extension,mask,Facial Mask,,,CoatPattern_locusE_mask.svg\n'
        '1,CoatPattern,ticking,no_ticking,,,,\n'
        '1,CoatPattern,brindle,not_brindle,,,,\n'
        '1,CoatPattern,merle,not_merle,,,,\n'
        '1,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '1,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '1,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '1,CoatType,Shedding Propensity,normal_shedding,Normal shedding,,,CoatType_Shedding_normal.svg\n'
        '1,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '1,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '1,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '2,CoatColor,black,black_leathers,Black (nose and paws),,#17181D,\n'
        '2,CoatColor,red,red,Red,,#7E341B,\n'
        '2,CoatPattern,agouti,tan_points_hidden,Tan Points (hidden),,,CoatPattern_locusA_tanPoints.svg\n'
        '2,CoatPattern,extension,rec_red,"No mask, grizzle, or domino patterns",,,\n'
        '2,CoatPattern,ticking,no_ticking,,,,\n'
        '2,CoatPattern,brindle,not_brindle,,,,\n'
        '2,CoatPattern,merle,not_merle,,,,\n'
        '2,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '2,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '2,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '2,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '2,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '2,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '2,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '3,CoatColor,black,black_solid,Black (solid coat),,#17181D,\n'
        '3,CoatColor,red,tan,Tan,,#A86B39,\n'
        '3,CoatPattern,agouti,tan_points_hidden,Tan Points (hidden),,,CoatPattern_locusA_tanPoints.svg\n'
        '3,CoatPattern,extension,mask_hidden,Facial Mask (hidden),,,CoatPattern_locusE_mask.svg\n'
        '3,CoatPattern,ticking,no_ticking,,,,\n'
        '3,CoatPattern,brindle,not_brindle,,,,\n'
        '3,CoatPattern,merle,not_merle,,,,\n'
        '3,CoatType,Coat Texture,wavy_coat,Wavy coat,,,CoatType_Curl_wavy.svg\n'
        '3,CoatType,Coat Length,long_coat,Long coat,,,CoatType_Length_long.svg\n'
        '3,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '3,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '3,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '3,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '3,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '4,CoatColor,black,black_solid,Black (solid coat),,#17181D,\n'
        '4,CoatColor,red,red,Red,,#7E341B,\n'
        '4,CoatPattern,agouti,sable_hidden,Sable (hidden),,,CoatPattern_locusA_sable.svg\n'
        '4,CoatPattern,extension,normal_extension_hidden,"No mask, grizzle, or domino patterns",,,\n'
        '4,CoatPattern,ticking,no_ticking,,,,\n'
        '4,CoatPattern,brindle,not_brindle,,,,\n'
        '4,CoatPattern,merle,not_merle,,,,\n'
        '4,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '4,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '4,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '4,CoatType,Shedding Propensity,normal_shedding,Normal shedding,,,CoatType_Shedding_normal.svg\n'
        '4,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '4,SpecialFeatures,Skeletal - Leg Length,short_legs_marker,Shortened leg length,,,SpecialFeatures_Limbs_short.svg\n'
        '4,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '5,CoatColor,black,black,Black,,#17181D,\n'
        '5,CoatColor,red,tan,Tan,,#A86B39,\n'
        '5,CoatPattern,agouti,tan_points,Tan Points,,,CoatPattern_locusA_tanPoints.svg\n'
        '5,CoatPattern,extension,mask,Facial Mask,,,CoatPattern_locusE_mask.svg\n'
        '5,CoatPattern,ticking,no_ticking,,,,\n'
        '5,CoatPattern,brindle,not_brindle,,,,\n'
        '5,CoatPattern,merle,not_merle,,,,\n'
        '5,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '5,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '5,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '5,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '5,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '5,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '5,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
    )

def genotype_table():
    return StringIO(
        # compared to legacy version, probabilities are always floats
        # '2,1,42300805,SpecialFeatures,leg_length,esr1,e3,G,*,1.0\n'
        # is now
        # '2,1,42300805,SpecialFeatures,leg_length,esr1,e3,T,G,1.0\n'
        # due to how multiple variants are resolved
        'sample,CHR,POS,tab,trait,locus,variant,A1,A2,CONF\n'
        '1,11,33317810,CoatColor,black,B,bc,T,T,1.0\n'
        '1,11,33326726,CoatColor,black,B,bd,ACCT,ACCT,1.0\n'
        '1,11,33326685,CoatColor,black,B,bs,C,C,1.0\n'
        '1,23,43969695,CoatColor,black,Co,co,G,G,1.0\n'
        '1,25,48121642,CoatColor,black,D,d1,G,G,1.0\n'
        '1,25,48150787,CoatColor,black,D,d2,G,G,1.0\n'
        '1,2,74746906,CoatColor,red,I,i1,T,A,1.0\n'
        '1,15,29840789,CoatColor,red,I,i2,C,C,1.0\n'
        '1,18,12910382,CoatColor,red,I,i3,T,T,1.0\n'
        '1,20,55850145,CoatColor,red,I,i4,C,C,1.0\n'
        '1,21,10864834,CoatColor,red,I,i5,G,G,1.0\n'
        '1,24,23393510,CoatPattern,agouti,A,Ay1,T,T,1.0\n'
        '1,24,23393514,CoatPattern,agouti,A,Ay2,A,A,1.0\n'
        '1,24,23365205,CoatPattern,agouti,A,at,C,C,1.0\n'
        '1,24,23393552,CoatPattern,agouti,A,a-,C,C,1.0\n'
        '1,24,23252763,CoatPattern,agouti,A,atm1,C,C,1.0\n'
        '1,24,23252753,CoatPattern,agouti,A,atm2,G,G,1.0\n'
        '1,5,63694460,CoatPattern,mask,E,EM,C,T,1.0\n'
        '1,5,63695017,CoatPattern,grizzle,E,Eg,C,C,1.0\n'
        '1,5,63694349,CoatPattern,northern_domino,E,EA,G,G,1.0\n'
        '1,5,63694334,CoatPattern,red,E,e1,G,A,1.0\n'
        '1,5,63695679,CoatPattern,red,E,e2,C,C,1.0\n'
        '1,16,58965448,CoatPattern,brindle,K,KB,TCCC,TCCC,1.0\n'
        '1,16,58999689,CoatPattern,brindle,K,Kb1,A,A,1.0\n'
        '1,16,58999678,CoatPattern,brindle,K,Kb2,GCTTCCCTAAAA,GCTTCCCTAAAA,1.0\n'
        '1,38,11165134,CoatPattern,roan,T,ta,G,G,1.0\n'
        '1,9,58530295,CoatPattern,merle,H,h,T,T,1.0\n'
        '1,27,2539211,CoatType,coat_texture,Cu,c1,C,C,1.0\n'
        '1,27,2543231,CoatType,coat_texture,Cu,c2,CTG,CTG,1.0\n'
        '1,32,4509367,CoatType,coat_length,L,Lh,G,G,1.0\n'
        '1,13,8491477,CoatType,furnishings,IC,ins,A,A,1.0\n'
        '1,1,24430748,CoatType,shedding,sp,sp,T,T,1.0\n'
        '1,28,24860187,CoatType,coat_layer,sc,1,T,T,1.0\n'
        '1,28,24870184,CoatType,coat_layer,sc,2,A,A,1.0\n'
        '1,10,48626862,SpecialFeatures,altitude,alt1,a1,G,G,1.0\n'
        '1,10,48630137,SpecialFeatures,altitude,alt2,a2,T,T,1.0\n'
        '1,10,48630153,SpecialFeatures,altitude,alt3,a3,G,G,1.0\n'
        '1,10,48633379,SpecialFeatures,altitude,alt4,a4,C,C,1.0\n'
        '1,18,44924848,SpecialFeatures,eye_color,alx4,blue,C,C,1.0\n'
        '1,18,20423056,SpecialFeatures,leg_length,fgf4c18gwas,snp,A,A,1.0\n'
        '1,1,42302476,SpecialFeatures,leg_length,esr1,e1,T,T,1.0\n'
        '1,1,42306135,SpecialFeatures,leg_length,esr1,e2,G,G,1.0\n'
        '1,1,42300805,SpecialFeatures,leg_length,esr1,e3,G,G,1.0\n'
        '1,1,54192143,SpecialFeatures,tail_length,BT,bt,NA,NA,1.0\n'
        '2,11,33317810,CoatColor,black,B,bc,T,T,1.0\n'
        '2,11,33326726,CoatColor,black,B,bd,ACCT,ACCT,1.0\n'
        '2,11,33326685,CoatColor,black,B,bs,C,C,1.0\n'
        '2,23,43969695,CoatColor,black,Co,co,G,G,1.0\n'
        '2,25,48121642,CoatColor,black,D,d1,G,G,1.0\n'
        '2,25,48150787,CoatColor,black,D,d2,G,G,1.0\n'
        '2,2,74746906,CoatColor,red,I,i1,A,A,1.0\n'
        '2,15,29840789,CoatColor,red,I,i2,T,C,1.0\n'
        '2,18,12910382,CoatColor,red,I,i3,T,T,1.0\n'
        '2,20,55850145,CoatColor,red,I,i4,C,C,1.0\n'
        '2,21,10864834,CoatColor,red,I,i5,G,A,1.0\n'
        '2,24,23393510,CoatPattern,agouti,A,Ay1,G,G,1.0\n'
        '2,24,23393514,CoatPattern,agouti,A,Ay2,G,G,1.0\n'
        '2,24,23365205,CoatPattern,agouti,A,at,T,T,1.0\n'
        '2,24,23393552,CoatPattern,agouti,A,a-,C,C,1.0\n'
        '2,24,23252763,CoatPattern,agouti,A,atm1,C,C,1.0\n'
        '2,24,23252753,CoatPattern,agouti,A,atm2,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,1.0\n'
        '2,5,63694460,CoatPattern,mask,E,EM,T,T,1.0\n'
        '2,5,63695017,CoatPattern,grizzle,E,Eg,C,C,1.0\n'
        '2,5,63694349,CoatPattern,northern_domino,E,EA,G,G,1.0\n'
        '2,5,63694334,CoatPattern,red,E,e1,A,A,1.0\n'
        '2,5,63695679,CoatPattern,red,E,e2,C,C,1.0\n'
        '2,16,58965448,CoatPattern,brindle,K,KB,T,T,1.0\n'
        '2,16,58999689,CoatPattern,brindle,K,Kb1,A,A,1.0\n'
        '2,16,58999678,CoatPattern,brindle,K,Kb2,GCTTCCCTAAAA,GCTTCCCTAAAA,1.0\n'
        '2,38,11165134,CoatPattern,roan,T,ta,G,G,1.0\n'
        '2,9,58530295,CoatPattern,merle,H,h,T,T,1.0\n'
        '2,27,2539211,CoatType,coat_texture,Cu,c1,C,C,1.0\n'
        '2,27,2543231,CoatType,coat_texture,Cu,c2,CTG,CTG,1.0\n'
        '2,32,4509367,CoatType,coat_length,L,Lh,G,G,1.0\n'
        '2,13,8491477,CoatType,furnishings,IC,ins,A,A,1.0\n'
        '2,1,24430748,CoatType,shedding,sp,sp,T,C,1.0\n'
        '2,28,24860187,CoatType,coat_layer,sc,1,T,T,1.0\n'
        '2,28,24870184,CoatType,coat_layer,sc,2,G,G,1.0\n'
        '2,10,48626862,SpecialFeatures,altitude,alt1,a1,G,G,1.0\n'
        '2,10,48630137,SpecialFeatures,altitude,alt2,a2,T,T,1.0\n'
        '2,10,48630153,SpecialFeatures,altitude,alt3,a3,G,G,1.0\n'
        '2,10,48633379,SpecialFeatures,altitude,alt4,a4,C,C,1.0\n'
        '2,18,44924848,SpecialFeatures,eye_color,alx4,blue,C,T,1.0\n'
        '2,18,20423056,SpecialFeatures,leg_length,fgf4c18gwas,snp,A,A,1.0\n'
        '2,1,42302476,SpecialFeatures,leg_length,esr1,e1,T,T,1.0\n'
        '2,1,42306135,SpecialFeatures,leg_length,esr1,e2,G,G,1.0\n'
        '2,1,42300805,SpecialFeatures,leg_length,esr1,e3,T,G,1.0\n'
        '2,1,54192143,SpecialFeatures,tail_length,BT,bt,G,G,1.0\n'
        '3,11,33317810,CoatColor,black,B,bc,T,T,1.0\n'
        '3,11,33326726,CoatColor,black,B,bd,ACCT,ACCT,1.0\n'
        '3,11,33326685,CoatColor,black,B,bs,C,C,1.0\n'
        '3,23,43969695,CoatColor,black,Co,co,G,G,1.0\n'
        '3,25,48121642,CoatColor,black,D,d1,G,G,1.0\n'
        '3,25,48150787,CoatColor,black,D,d2,G,G,1.0\n'
        '3,2,74746906,CoatColor,red,I,i1,T,A,1.0\n'
        '3,15,29840789,CoatColor,red,I,i2,T,C,1.0\n'
        '3,18,12910382,CoatColor,red,I,i3,T,T,1.0\n'
        '3,20,55850145,CoatColor,red,I,i4,C,C,1.0\n'
        '3,21,10864834,CoatColor,red,I,i5,G,A,1.0\n'
        '3,24,23393510,CoatPattern,agouti,A,Ay1,G,G,1.0\n'
        '3,24,23393514,CoatPattern,agouti,A,Ay2,G,G,1.0\n'
        '3,24,23365205,CoatPattern,agouti,A,at,T,T,1.0\n'
        '3,24,23393552,CoatPattern,agouti,A,a-,C,T,1.0\n'
        '3,24,23252763,CoatPattern,agouti,A,atm1,C,C,1.0\n'
        '3,24,23252753,CoatPattern,agouti,A,atm2,G,G,1.0\n'
        '3,5,63694460,CoatPattern,mask,E,EM,C,C,1.0\n'
        '3,5,63695017,CoatPattern,grizzle,E,Eg,C,C,1.0\n'
        '3,5,63694349,CoatPattern,northern_domino,E,EA,G,G,1.0\n'
        '3,5,63694334,CoatPattern,red,E,e1,G,G,1.0\n'
        '3,5,63695679,CoatPattern,red,E,e2,C,C,1.0\n'
        '3,16,58965448,CoatPattern,brindle,K,KB,TCCC,T,1.0\n'
        '3,16,58999689,CoatPattern,brindle,K,Kb1,A,A,1.0\n'
        '3,16,58999678,CoatPattern,brindle,K,Kb2,GCTTCCCTAAAA,GCTTCCCTAAAA,1.0\n'
        '3,38,11165134,CoatPattern,roan,T,ta,G,A,1.0\n'
        '3,9,58530295,CoatPattern,merle,H,h,T,T,1.0\n'
        '3,27,2539211,CoatType,coat_texture,Cu,c1,C,T,1.0\n'
        '3,27,2543231,CoatType,coat_texture,Cu,c2,CTG,CTG,1.0\n'
        '3,32,4509367,CoatType,coat_length,L,Lh,T,T,1.0\n'
        '3,13,8491477,CoatType,furnishings,IC,ins,A,A,1.0\n'
        '3,1,24430748,CoatType,shedding,sp,sp,C,C,1.0\n'
        '3,28,24860187,CoatType,coat_layer,sc,1,C,C,1.0\n'
        '3,28,24870184,CoatType,coat_layer,sc,2,G,G,1.0\n'
        '3,10,48626862,SpecialFeatures,altitude,alt1,a1,G,G,1.0\n'
        '3,10,48630137,SpecialFeatures,altitude,alt2,a2,T,T,1.0\n'
        '3,10,48630153,SpecialFeatures,altitude,alt3,a3,G,G,1.0\n'
        '3,10,48633379,SpecialFeatures,altitude,alt4,a4,C,C,1.0\n'
        '3,18,44924848,SpecialFeatures,eye_color,alx4,blue,C,C,1.0\n'
        '3,18,20423056,SpecialFeatures,leg_length,fgf4c18gwas,snp,A,A,1.0\n'
        '3,1,42302476,SpecialFeatures,leg_length,esr1,e1,T,T,1.0\n'
        '3,1,42306135,SpecialFeatures,leg_length,esr1,e2,G,G,1.0\n'
        '3,1,42300805,SpecialFeatures,leg_length,esr1,e3,G,G,1.0\n'
        '3,1,54192143,SpecialFeatures,tail_length,BT,bt,G,G,1.0\n'
        '4,11,33317810,CoatColor,black,B,bc,T,T,1.0\n'
        '4,11,33326726,CoatColor,black,B,bd,ACCT,ACCT,1.0\n'
        '4,11,33326685,CoatColor,black,B,bs,C,T,1.0\n'
        '4,23,43969695,CoatColor,black,Co,co,G,G,1.0\n'
        '4,25,48121642,CoatColor,black,D,d1,G,G,1.0\n'
        '4,25,48150787,CoatColor,black,D,d2,G,G,1.0\n'
        '4,2,74746906,CoatColor,red,I,i1,A,A,1.0\n'
        '4,15,29840789,CoatColor,red,I,i2,T,C,1.0\n'
        '4,18,12910382,CoatColor,red,I,i3,T,T,1.0\n'
        '4,20,55850145,CoatColor,red,I,i4,C,C,1.0\n'
        '4,21,10864834,CoatColor,red,I,i5,A,A,1.0\n'
        '4,24,23393510,CoatPattern,agouti,A,Ay1,T,G,1.0\n'
        '4,24,23393514,CoatPattern,agouti,A,Ay2,A,G,1.0\n'
        '4,24,23365205,CoatPattern,agouti,A,at,C,T,1.0\n'
        '4,24,23393552,CoatPattern,agouti,A,a-,C,C,1.0\n'
        '4,24,23252763,CoatPattern,agouti,A,atm1,C,C,1.0\n'
        '4,24,23252753,CoatPattern,agouti,A,atm2,G,GTCCCCAGGTCAGAGTT,1.0\n'
        '4,5,63694460,CoatPattern,mask,E,EM,T,T,1.0\n'
        '4,5,63695017,CoatPattern,grizzle,E,Eg,C,C,1.0\n'
        '4,5,63694349,CoatPattern,northern_domino,E,EA,G,G,1.0\n'
        '4,5,63694334,CoatPattern,red,E,e1,G,A,1.0\n'
        '4,5,63695679,CoatPattern,red,E,e2,C,C,1.0\n'
        '4,16,58965448,CoatPattern,brindle,K,KB,TCCC,T,1.0\n'
        '4,16,58999689,CoatPattern,brindle,K,Kb1,A,A,1.0\n'
        '4,16,58999678,CoatPattern,brindle,K,Kb2,GCTTCCCTAAAA,GCTTCCCTAAAA,1.0\n'
        '4,38,11165134,CoatPattern,roan,T,ta,G,G,1.0\n'
        '4,9,58530295,CoatPattern,merle,H,h,T,T,1.0\n'
        '4,27,2539211,CoatType,coat_texture,Cu,c1,C,C,1.0\n'
        '4,27,2543231,CoatType,coat_texture,Cu,c2,CTG,CTG,1.0\n'
        '4,32,4509367,CoatType,coat_length,L,Lh,G,G,1.0\n'
        '4,13,8491477,CoatType,furnishings,IC,ins,A,A,1.0\n'
        '4,1,24430748,CoatType,shedding,sp,sp,T,T,1.0\n'
        '4,28,24860187,CoatType,coat_layer,sc,1,T,T,1.0\n'
        '4,28,24870184,CoatType,coat_layer,sc,2,A,A,1.0\n'
        '4,10,48626862,SpecialFeatures,altitude,alt1,a1,G,G,1.0\n'
        '4,10,48630137,SpecialFeatures,altitude,alt2,a2,T,T,1.0\n'
        '4,10,48630153,SpecialFeatures,altitude,alt3,a3,G,G,1.0\n'
        '4,10,48633379,SpecialFeatures,altitude,alt4,a4,C,C,1.0\n'
        '4,18,44924848,SpecialFeatures,eye_color,alx4,blue,C,C,1.0\n'
        '4,18,20423056,SpecialFeatures,leg_length,fgf4c18gwas,snp,G,G,1.0\n'
        '4,1,42302476,SpecialFeatures,leg_length,esr1,e1,T,T,1.0\n'
        '4,1,42306135,SpecialFeatures,leg_length,esr1,e2,G,G,1.0\n'
        '4,1,42300805,SpecialFeatures,leg_length,esr1,e3,G,G,1.0\n'
        '4,1,54192143,SpecialFeatures,tail_length,BT,bt,G,G,1.0\n'
        '5,11,33317810,CoatColor,black,B,bc,T,T,1.0\n'
        '5,11,33326726,CoatColor,black,B,bd,ACCT,ACCT,1.0\n'
        '5,11,33326685,CoatColor,black,B,bs,C,C,1.0\n'
        '5,23,43969695,CoatColor,black,Co,co,G,G,1.0\n'
        '5,25,48121642,CoatColor,black,D,d1,G,G,1.0\n'
        '5,25,48150787,CoatColor,black,D,d2,G,G,1.0\n'
        '5,2,74746906,CoatColor,red,I,i1,T,T,1.0\n'
        '5,15,29840789,CoatColor,red,I,i2,C,C,1.0\n'
        '5,18,12910382,CoatColor,red,I,i3,T,C,1.0\n'
        '5,20,55850145,CoatColor,red,I,i4,C,C,1.0\n'
        '5,21,10864834,CoatColor,red,I,i5,G,G,1.0\n'
        '5,24,23393510,CoatPattern,agouti,A,Ay1,G,G,1.0\n'
        '5,24,23393514,CoatPattern,agouti,A,Ay2,G,G,1.0\n'
        '5,24,23365205,CoatPattern,agouti,A,at,T,T,1.0\n'
        '5,24,23393552,CoatPattern,agouti,A,a-,C,C,1.0\n'
        '5,24,23252763,CoatPattern,agouti,A,atm1,C,C,1.0\n'
        '5,24,23252753,CoatPattern,agouti,A,atm2,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,1.0\n'
        '5,5,63694460,CoatPattern,mask,E,EM,C,T,1.0\n'
        '5,5,63695017,CoatPattern,grizzle,E,Eg,C,C,1.0\n'
        '5,5,63694349,CoatPattern,northern_domino,E,EA,G,G,1.0\n'
        '5,5,63694334,CoatPattern,red,E,e1,G,G,1.0\n'
        '5,5,63695679,CoatPattern,red,E,e2,C,C,1.0\n'
        '5,16,58965448,CoatPattern,brindle,K,KB,TCCC,TCCC,1.0\n'
        '5,16,58999689,CoatPattern,brindle,K,Kb1,A,A,1.0\n'
        '5,16,58999678,CoatPattern,brindle,K,Kb2,GCTTCCCTAAAA,GCTTCCCTAAAA,1.0\n'
        '5,38,11165134,CoatPattern,roan,T,ta,G,A,1.0\n'
        '5,9,58530295,CoatPattern,merle,H,h,T,T,1.0\n'
        '5,27,2539211,CoatType,coat_texture,Cu,c1,C,C,1.0\n'
        '5,27,2543231,CoatType,coat_texture,Cu,c2,CTG,CTG,1.0\n'
        '5,32,4509367,CoatType,coat_length,L,Lh,G,G,1.0\n'
        '5,13,8491477,CoatType,furnishings,IC,ins,A,A,1.0\n'
        '5,1,24430748,CoatType,shedding,sp,sp,T,C,1.0\n'
        '5,28,24860187,CoatType,coat_layer,sc,1,T,C,1.0\n'
        '5,28,24870184,CoatType,coat_layer,sc,2,A,G,1.0\n'
        '5,10,48626862,SpecialFeatures,altitude,alt1,a1,G,G,1.0\n'
        '5,10,48630137,SpecialFeatures,altitude,alt2,a2,T,T,1.0\n'
        '5,10,48630153,SpecialFeatures,altitude,alt3,a3,G,G,1.0\n'
        '5,10,48633379,SpecialFeatures,altitude,alt4,a4,C,C,1.0\n'
        '5,18,44924848,SpecialFeatures,eye_color,alx4,blue,C,C,1.0\n'
        '5,18,20423056,SpecialFeatures,leg_length,fgf4c18gwas,snp,A,A,1.0\n'
        '5,1,42302476,SpecialFeatures,leg_length,esr1,e1,T,T,1.0\n'
        '5,1,42306135,SpecialFeatures,leg_length,esr1,e2,G,G,1.0\n'
        '5,1,42300805,SpecialFeatures,leg_length,esr1,e3,G,G,1.0\n'
        '5,1,54192143,SpecialFeatures,tail_length,BT,bt,G,G,1.0'
    )

def phenotype_table():
    return StringIO(
        'sample,tab,trait,result,string,value,color,image\n'
        '1,CoatColor,black,black,Black,,#17181D,\n'
        '1,CoatColor,red,tan,Tan,,#A86B39,\n'
        '1,CoatPattern,agouti,sable,Sable,,,CoatPattern_locusA_sable.svg\n'
        '1,CoatPattern,extension,mask,Facial Mask,,,CoatPattern_locusE_mask.svg\n'
        '1,CoatPattern,ticking,no_ticking,,,,\n'
        '1,CoatPattern,brindle,not_brindle,,,,\n'
        '1,CoatPattern,merle,not_merle,,,,\n'
        '1,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '1,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '1,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '1,CoatType,Shedding Propensity,normal_shedding,Normal shedding,,,CoatType_Shedding_normal.svg\n'
        '1,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '1,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '1,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '2,CoatColor,black,black_leathers,Black (nose and paws),,#17181D,\n'
        '2,CoatColor,red,red,Red,,#7E341B,\n'
        '2,CoatPattern,agouti,tan_points_hidden,Tan Points (hidden),,,CoatPattern_locusA_tanPoints.svg\n'
        '2,CoatPattern,extension,rec_red,"No mask, grizzle, or domino patterns",,,\n'
        '2,CoatPattern,ticking,no_ticking,,,,\n'
        '2,CoatPattern,brindle,not_brindle,,,,\n'
        '2,CoatPattern,merle,not_merle,,,,\n'
        '2,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '2,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '2,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '2,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '2,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '2,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '2,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '3,CoatColor,black,black_solid,Black (solid coat),,#17181D,\n'
        '3,CoatColor,red,tan,Tan,,#A86B39,\n'
        '3,CoatPattern,agouti,tan_points_hidden,Tan Points (hidden),,,CoatPattern_locusA_tanPoints.svg\n'
        '3,CoatPattern,extension,mask_hidden,Facial Mask (hidden),,,CoatPattern_locusE_mask.svg\n'
        '3,CoatPattern,ticking,no_ticking,,,,\n'
        '3,CoatPattern,brindle,not_brindle,,,,\n'
        '3,CoatPattern,merle,not_merle,,,,\n'
        '3,CoatType,Coat Texture,wavy_coat,Wavy coat,,,CoatType_Curl_wavy.svg\n'
        '3,CoatType,Coat Length,long_coat,Long coat,,,CoatType_Length_long.svg\n'
        '3,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '3,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '3,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '3,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '3,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '4,CoatColor,black,black_solid,Black (solid coat),,#17181D,\n'
        '4,CoatColor,red,red,Red,,#7E341B,\n'
        '4,CoatPattern,agouti,sable_hidden,Sable (hidden),,,CoatPattern_locusA_sable.svg\n'
        '4,CoatPattern,extension,normal_extension_hidden,"No mask, grizzle, or domino patterns",,,\n'
        '4,CoatPattern,ticking,no_ticking,,,,\n'
        '4,CoatPattern,brindle,not_brindle,,,,\n'
        '4,CoatPattern,merle,not_merle,,,,\n'
        '4,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '4,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '4,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '4,CoatType,Shedding Propensity,normal_shedding,Normal shedding,,,CoatType_Shedding_normal.svg\n'
        '4,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '4,SpecialFeatures,Skeletal - Leg Length,short_legs_marker,Shortened leg length,,,SpecialFeatures_Limbs_short.svg\n'
        '4,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '5,CoatColor,black,black,Black,,#17181D,\n'
        '5,CoatColor,red,tan,Tan,,#A86B39,\n'
        '5,CoatPattern,agouti,tan_points,Tan Points,,,CoatPattern_locusA_tanPoints.svg\n'
        '5,CoatPattern,extension,mask,Facial Mask,,,CoatPattern_locusE_mask.svg\n'
        '5,CoatPattern,ticking,no_ticking,,,,\n'
        '5,CoatPattern,brindle,not_brindle,,,,\n'
        '5,CoatPattern,merle,not_merle,,,,\n'
        '5,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '5,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '5,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '5,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '5,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '5,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '5,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg'
    )

def json_table():
    return StringIO(
        'sample,tab,name,gene,firstCopy,secondCopy,possibleAlleles,effect\n'
        '1,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '2,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '3,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '4,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '5,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '1,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '2,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '3,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '4,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '5,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '1,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '2,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '3,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '4,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,T,C & T,0\n'
        '5,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '1,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '2,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '3,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '4,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '5,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '1,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '2,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '3,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '4,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '5,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '1,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '2,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '3,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '4,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '5,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '1,CoatColor,Red intensity - marker 1,lincRNA,T,A,T & A,0\n'
        '2,CoatColor,Red intensity - marker 1,lincRNA,A,A,T & A,0\n'
        '3,CoatColor,Red intensity - marker 1,lincRNA,T,A,T & A,0\n'
        '4,CoatColor,Red intensity - marker 1,lincRNA,A,A,T & A,0\n'
        '5,CoatColor,Red intensity - marker 1,lincRNA,T,T,T & A,0\n'
        '1,CoatColor,Red intensity - marker 2,intergenic,C,C,T & C,0\n'
        '2,CoatColor,Red intensity - marker 2,intergenic,T,C,T & C,0\n'
        '3,CoatColor,Red intensity - marker 2,intergenic,T,C,T & C,0\n'
        '4,CoatColor,Red intensity - marker 2,intergenic,T,C,T & C,0\n'
        '5,CoatColor,Red intensity - marker 2,intergenic,C,C,T & C,0\n'
        '1,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '2,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '3,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '4,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '5,CoatColor,Red intensity - marker 3,SLC264A,T,C,T & C,0\n'
        '1,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '2,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '3,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '4,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '5,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '1,CoatColor,Red intensity - marker 5,TYR,G,G,G & A,0\n'
        '2,CoatColor,Red intensity - marker 5,TYR,G,A,G & A,0\n'
        '3,CoatColor,Red intensity - marker 5,TYR,G,A,G & A,0\n'
        '4,CoatColor,Red intensity - marker 5,TYR,A,A,G & A,0\n'
        '5,CoatColor,Red intensity - marker 5,TYR,G,G,G & A,0\n'
        '1,CoatPattern,Sable - variant p.(A82S),ASIP,T,T,G & T,0\n'
        '2,CoatPattern,Sable - variant p.(A82S),ASIP,G,G,G & T,0\n'
        '3,CoatPattern,Sable - variant p.(A82S),ASIP,G,G,G & T,0\n'
        '4,CoatPattern,Sable - variant p.(A82S),ASIP,T,G,G & T,0\n'
        '5,CoatPattern,Sable - variant p.(A82S),ASIP,G,G,G & T,0\n'
        '1,CoatPattern,Sable - variant p.(R83H),ASIP,A,A,G & A,0\n'
        '2,CoatPattern,Sable - variant p.(R83H),ASIP,G,G,G & A,0\n'
        '3,CoatPattern,Sable - variant p.(R83H),ASIP,G,G,G & A,0\n'
        '4,CoatPattern,Sable - variant p.(R83H),ASIP,A,G,G & A,0\n'
        '5,CoatPattern,Sable - variant p.(R83H),ASIP,G,G,G & A,0\n'
        '1,CoatPattern,Tan points - marker,ASIP,C,C,C & T,0\n'
        '2,CoatPattern,Tan points - marker,ASIP,T,T,C & T,0\n'
        '3,CoatPattern,Tan points - marker,ASIP,T,T,C & T,0\n'
        '4,CoatPattern,Tan points - marker,ASIP,C,T,C & T,0\n'
        '5,CoatPattern,Tan points - marker,ASIP,T,T,C & T,0\n'
        '1,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,C,C & T,0\n'
        '2,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,C,C & T,0\n'
        '3,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,T,C & T,0\n'
        '4,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,C,C & T,0\n'
        '5,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,C,C & T,0\n'
        '1,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '2,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '3,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '4,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '5,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '1,CoatPattern,Saddle - marker 2,RALY,G,G,GTCCCCAGGTCAGAGTT & G,0\n'
        '2,CoatPattern,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT & G,0\n'
        '3,CoatPattern,Saddle - marker 2,RALY,G,G,GTCCCCAGGTCAGAGTT & G,0\n'
        '4,CoatPattern,Saddle - marker 2,RALY,G,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT & G,0\n'
        '5,CoatPattern,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT & G,0\n'
        '1,CoatPattern,Facial mask - variant p.(M264V),MC1R,C,T,T & C,0\n'
        '2,CoatPattern,Facial mask - variant p.(M264V),MC1R,T,T,T & C,0\n'
        '3,CoatPattern,Facial mask - variant p.(M264V),MC1R,C,C,T & C,0\n'
        '4,CoatPattern,Facial mask - variant p.(M264V),MC1R,T,T,T & C,0\n'
        '5,CoatPattern,Facial mask - variant p.(M264V),MC1R,C,T,T & C,0\n'
        '1,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '2,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '3,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '4,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '5,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '1,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '2,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '3,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '4,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '5,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '1,CoatPattern,Recessive red - variant p.(R306*),MC1R,G,A,G & A,0\n'
        '2,CoatPattern,Recessive red - variant p.(R306*),MC1R,A,A,G & A,0\n'
        '3,CoatPattern,Recessive red - variant p.(R306*),MC1R,G,G,G & A,0\n'
        '4,CoatPattern,Recessive red - variant p.(R306*),MC1R,G,A,G & A,0\n'
        '5,CoatPattern,Recessive red - variant p.(R306*),MC1R,G,G,G & A,0\n'
        '1,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '2,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '3,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '4,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '5,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '1,CoatPattern,Dominant black - variant p.(G23del),CBD103,TCCC,TCCC,TCCC & T,0\n'
        '2,CoatPattern,Dominant black - variant p.(G23del),CBD103,T,T,TCCC & T,0\n'
        '3,CoatPattern,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC & T,0\n'
        '4,CoatPattern,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC & T,0\n'
        '5,CoatPattern,Dominant black - variant p.(G23del),CBD103,TCCC,TCCC,TCCC & T,0\n'
        '1,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '2,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '3,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '4,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '5,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '1,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '2,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '3,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '4,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '5,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '1,CoatPattern,Ticking - marker,USH2A,G,G,G & A,0\n'
        '2,CoatPattern,Ticking - marker,USH2A,G,G,G & A,0\n'
        '3,CoatPattern,Ticking - marker,USH2A,G,A,G & A,0\n'
        '4,CoatPattern,Ticking - marker,USH2A,G,G,G & A,0\n'
        '5,CoatPattern,Ticking - marker,USH2A,G,A,G & A,0\n'
        '1,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '2,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '3,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '4,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '5,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '1,CoatType,Curly coat - variant p.(R151W),KRT71,C,C,C & T,0\n'
        '2,CoatType,Curly coat - variant p.(R151W),KRT71,C,C,C & T,0\n'
        '3,CoatType,Curly coat - variant p.(R151W),KRT71,C,T,C & T,0\n'
        '4,CoatType,Curly coat - variant p.(R151W),KRT71,C,C,C & T,0\n'
        '5,CoatType,Curly coat - variant p.(R151W),KRT71,C,C,C & T,0\n'
        '1,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '2,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '3,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '4,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '5,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '1,CoatType,Long coat - variant p.(C95F),FGF5,G,G,G & T,0\n'
        '2,CoatType,Long coat - variant p.(C95F),FGF5,G,G,G & T,0\n'
        '3,CoatType,Long coat - variant p.(C95F),FGF5,T,T,G & T,0\n'
        '4,CoatType,Long coat - variant p.(C95F),FGF5,G,G,G & T,0\n'
        '5,CoatType,Long coat - variant p.(C95F),FGF5,G,G,G & T,0\n'
        '1,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '2,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '3,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '4,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '5,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '1,CoatType,Shedding propensity - variant p.(A237T),MC5R,T,T,T & C,0\n'
        '2,CoatType,Shedding propensity - variant p.(A237T),MC5R,T,C,T & C,0\n'
        '3,CoatType,Shedding propensity - variant p.(A237T),MC5R,C,C,T & C,0\n'
        '4,CoatType,Shedding propensity - variant p.(A237T),MC5R,T,T,T & C,0\n'
        '5,CoatType,Shedding propensity - variant p.(A237T),MC5R,T,C,T & C,0\n'
        '1,CoatType,Single-layer coat - marker 1,ADRB1-AU1,T,T,C & T,0\n'
        '2,CoatType,Single-layer coat - marker 1,ADRB1-AU1,T,T,C & T,0\n'
        '3,CoatType,Single-layer coat - marker 1,ADRB1-AU1,C,C,C & T,0\n'
        '4,CoatType,Single-layer coat - marker 1,ADRB1-AU1,T,T,C & T,0\n'
        '5,CoatType,Single-layer coat - marker 1,ADRB1-AU1,T,C,C & T,0\n'
        '1,CoatType,Single-layer coat - marker 2,ADRB1-AU1,A,A,G & A,0\n'
        '2,CoatType,Single-layer coat - marker 2,ADRB1-AU1,G,G,G & A,0\n'
        '3,CoatType,Single-layer coat - marker 2,ADRB1-AU1,G,G,G & A,0\n'
        '4,CoatType,Single-layer coat - marker 2,ADRB1-AU1,A,A,G & A,0\n'
        '5,CoatType,Single-layer coat - marker 2,ADRB1-AU1,A,G,G & A,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '2,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '3,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '4,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '5,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '2,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '3,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '4,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '5,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '2,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '3,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '4,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '5,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '2,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '3,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '4,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '5,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '1,SpecialFeatures,Blue eyes - marker,ALX4,C,C,C & T,0\n'
        '2,SpecialFeatures,Blue eyes - marker,ALX4,C,T,C & T,0\n'
        '3,SpecialFeatures,Blue eyes - marker,ALX4,C,C,C & T,0\n'
        '4,SpecialFeatures,Blue eyes - marker,ALX4,C,C,C & T,0\n'
        '5,SpecialFeatures,Blue eyes - marker,ALX4,C,C,C & T,0\n'
        '1,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '2,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '3,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '4,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,G,G,A & G,0\n'
        '5,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '1,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '2,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '3,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '4,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '5,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '1,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '2,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '3,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '4,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '5,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '1,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '2,SpecialFeatures,Long legs - marker 3,ESR1,T,G,T & G,0\n'
        '3,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '4,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '5,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '2,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
        '3,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
        '4,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
        '5,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0'
    )

def trail_geno():
    return StringIO(
        'Sample ID,tab,trait,Trait,Gene,Normal Version,Variant Version,First Copy,Second Copy,Description\n'
        '1,CoatColor,black,Liver - variant p.(C41S),TYRP1,T,A,T,T,Alters eumelanin pigment to a brown hue\n'
        '2,CoatColor,black,Liver - variant p.(C41S),TYRP1,T,A,T,T,Alters eumelanin pigment to a brown hue\n'
        '3,CoatColor,black,Liver - variant p.(C41S),TYRP1,T,A,T,T,Alters eumelanin pigment to a brown hue\n'
        '4,CoatColor,black,Liver - variant p.(C41S),TYRP1,T,A,T,T,Alters eumelanin pigment to a brown hue\n'
        '5,CoatColor,black,Liver - variant p.(C41S),TYRP1,T,A,T,T,Alters eumelanin pigment to a brown hue\n'
        '1,CoatColor,black,Liver - variant p.(P345del),TYRP1,ACCT,A,ACCT,ACCT,Alters eumelanin pigment to a brown hue\n'
        '2,CoatColor,black,Liver - variant p.(P345del),TYRP1,ACCT,A,ACCT,ACCT,Alters eumelanin pigment to a brown hue\n'
        '3,CoatColor,black,Liver - variant p.(P345del),TYRP1,ACCT,A,ACCT,ACCT,Alters eumelanin pigment to a brown hue\n'
        '4,CoatColor,black,Liver - variant p.(P345del),TYRP1,ACCT,A,ACCT,ACCT,Alters eumelanin pigment to a brown hue\n'
        '5,CoatColor,black,Liver - variant p.(P345del),TYRP1,ACCT,A,ACCT,ACCT,Alters eumelanin pigment to a brown hue\n'
        '1,CoatColor,black,Liver - variant p.(Gln331*),TYRP1,C,T,C,C,Alters eumelanin pigment to a brown hue\n'
        '2,CoatColor,black,Liver - variant p.(Gln331*),TYRP1,C,T,C,C,Alters eumelanin pigment to a brown hue\n'
        '3,CoatColor,black,Liver - variant p.(Gln331*),TYRP1,C,T,C,C,Alters eumelanin pigment to a brown hue\n'
        '4,CoatColor,black,Liver - variant p.(Gln331*),TYRP1,C,T,C,T,Alters eumelanin pigment to a brown hue\n'
        '5,CoatColor,black,Liver - variant p.(Gln331*),TYRP1,C,T,C,C,Alters eumelanin pigment to a brown hue\n'
        '1,CoatColor,black,Cocoa - variant p.(T807*),HPS3,G,A,G,G,Alters eumelanin pigment to a deep brown hue\n'
        '2,CoatColor,black,Cocoa - variant p.(T807*),HPS3,G,A,G,G,Alters eumelanin pigment to a deep brown hue\n'
        '3,CoatColor,black,Cocoa - variant p.(T807*),HPS3,G,A,G,G,Alters eumelanin pigment to a deep brown hue\n'
        '4,CoatColor,black,Cocoa - variant p.(T807*),HPS3,G,A,G,G,Alters eumelanin pigment to a deep brown hue\n'
        '5,CoatColor,black,Cocoa - variant p.(T807*),HPS3,G,A,G,G,Alters eumelanin pigment to a deep brown hue\n'
        '1,CoatColor,black,Dilution - splice variant,MLPH,G,A,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '2,CoatColor,black,Dilution - splice variant,MLPH,G,A,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '3,CoatColor,black,Dilution - splice variant,MLPH,G,A,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '4,CoatColor,black,Dilution - splice variant,MLPH,G,A,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '5,CoatColor,black,Dilution - splice variant,MLPH,G,A,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '1,CoatColor,black,Dilution - variant p.(Q235H),MLPH,G,C,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '2,CoatColor,black,Dilution - variant p.(Q235H),MLPH,G,C,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '3,CoatColor,black,Dilution - variant p.(Q235H),MLPH,G,C,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '4,CoatColor,black,Dilution - variant p.(Q235H),MLPH,G,C,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '5,CoatColor,black,Dilution - variant p.(Q235H),MLPH,G,C,G,G,Dilutes eumelanin pigment to a lighter hue\n'
        '1,CoatColor,red,Red intensity - marker 1,lincRNA,T,A,T,A,Influences deepness of pheomelanin pigment\n'
        '2,CoatColor,red,Red intensity - marker 1,lincRNA,T,A,A,A,Influences deepness of pheomelanin pigment\n'
        '3,CoatColor,red,Red intensity - marker 1,lincRNA,T,A,T,A,Influences deepness of pheomelanin pigment\n'
        '4,CoatColor,red,Red intensity - marker 1,lincRNA,T,A,A,A,Influences deepness of pheomelanin pigment\n'
        '5,CoatColor,red,Red intensity - marker 1,lincRNA,T,A,T,T,Influences deepness of pheomelanin pigment\n'
        '1,CoatColor,red,Red intensity - marker 2,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '2,CoatColor,red,Red intensity - marker 2,intergenic,T,C,T,C,Influences deepness of pheomelanin pigment\n'
        '3,CoatColor,red,Red intensity - marker 2,intergenic,T,C,T,C,Influences deepness of pheomelanin pigment\n'
        '4,CoatColor,red,Red intensity - marker 2,intergenic,T,C,T,C,Influences deepness of pheomelanin pigment\n'
        '5,CoatColor,red,Red intensity - marker 2,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '1,CoatColor,red,Red intensity - marker 3,SLC264A,T,C,T,T,Influences deepness of pheomelanin pigment\n'
        '2,CoatColor,red,Red intensity - marker 3,SLC264A,T,C,T,T,Influences deepness of pheomelanin pigment\n'
        '3,CoatColor,red,Red intensity - marker 3,SLC264A,T,C,T,T,Influences deepness of pheomelanin pigment\n'
        '4,CoatColor,red,Red intensity - marker 3,SLC264A,T,C,T,T,Influences deepness of pheomelanin pigment\n'
        '5,CoatColor,red,Red intensity - marker 3,SLC264A,T,C,T,C,Influences deepness of pheomelanin pigment\n'
        '1,CoatColor,red,Red intensity - marker 4,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '2,CoatColor,red,Red intensity - marker 4,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '3,CoatColor,red,Red intensity - marker 4,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '4,CoatColor,red,Red intensity - marker 4,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '5,CoatColor,red,Red intensity - marker 4,intergenic,T,C,C,C,Influences deepness of pheomelanin pigment\n'
        '1,CoatColor,red,Red intensity - marker 5,TYR,G,A,G,G,Influences deepness of pheomelanin pigment\n'
        '2,CoatColor,red,Red intensity - marker 5,TYR,G,A,G,A,Influences deepness of pheomelanin pigment\n'
        '3,CoatColor,red,Red intensity - marker 5,TYR,G,A,G,A,Influences deepness of pheomelanin pigment\n'
        '4,CoatColor,red,Red intensity - marker 5,TYR,G,A,A,A,Influences deepness of pheomelanin pigment\n'
        '5,CoatColor,red,Red intensity - marker 5,TYR,G,A,G,G,Influences deepness of pheomelanin pigment\n'
        '1,CoatPattern,agouti,Sable - variant p.(A82S),ASIP,G,T,T,T,Determines if a dog has sable pattern\n'
        '2,CoatPattern,agouti,Sable - variant p.(A82S),ASIP,G,T,G,G,Determines if a dog has sable pattern\n'
        '3,CoatPattern,agouti,Sable - variant p.(A82S),ASIP,G,T,G,G,Determines if a dog has sable pattern\n'
        '4,CoatPattern,agouti,Sable - variant p.(A82S),ASIP,G,T,T,G,Determines if a dog has sable pattern\n'
        '5,CoatPattern,agouti,Sable - variant p.(A82S),ASIP,G,T,G,G,Determines if a dog has sable pattern\n'
        '1,CoatPattern,agouti,Sable - variant p.(R83H),ASIP,G,A,A,A,Determines if a dog has sable pattern\n'
        '2,CoatPattern,agouti,Sable - variant p.(R83H),ASIP,G,A,G,G,Determines if a dog has sable pattern\n'
        '3,CoatPattern,agouti,Sable - variant p.(R83H),ASIP,G,A,G,G,Determines if a dog has sable pattern\n'
        '4,CoatPattern,agouti,Sable - variant p.(R83H),ASIP,G,A,A,G,Determines if a dog has sable pattern\n'
        '5,CoatPattern,agouti,Sable - variant p.(R83H),ASIP,G,A,G,G,Determines if a dog has sable pattern\n'
        '1,CoatPattern,agouti,Tan points - marker,ASIP,C,T,C,C,Determines if a dog has tan point pattern\n'
        '2,CoatPattern,agouti,Tan points - marker,ASIP,C,T,T,T,Determines if a dog has tan point pattern\n'
        '3,CoatPattern,agouti,Tan points - marker,ASIP,C,T,T,T,Determines if a dog has tan point pattern\n'
        '4,CoatPattern,agouti,Tan points - marker,ASIP,C,T,C,T,Determines if a dog has tan point pattern\n'
        '5,CoatPattern,agouti,Tan points - marker,ASIP,C,T,T,T,Determines if a dog has tan point pattern\n'
        '1,CoatPattern,agouti,Recessive black - variant p.(R96C),ASIP,C,T,C,C,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '2,CoatPattern,agouti,Recessive black - variant p.(R96C),ASIP,C,T,C,C,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '3,CoatPattern,agouti,Recessive black - variant p.(R96C),ASIP,C,T,C,T,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '4,CoatPattern,agouti,Recessive black - variant p.(R96C),ASIP,C,T,C,C,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '5,CoatPattern,agouti,Recessive black - variant p.(R96C),ASIP,C,T,C,C,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '1,CoatPattern,agouti,Saddle - marker 1,RALY,CAGAGTTTCCCCAGGT,C,C,C,Allows for saddle pattern to supercede tan points\n'
        '2,CoatPattern,agouti,Saddle - marker 1,RALY,CAGAGTTTCCCCAGGT,C,C,C,Allows for saddle pattern to supercede tan points\n'
        '3,CoatPattern,agouti,Saddle - marker 1,RALY,CAGAGTTTCCCCAGGT,C,C,C,Allows for saddle pattern to supercede tan points\n'
        '4,CoatPattern,agouti,Saddle - marker 1,RALY,CAGAGTTTCCCCAGGT,C,C,C,Allows for saddle pattern to supercede tan points\n'
        '5,CoatPattern,agouti,Saddle - marker 1,RALY,CAGAGTTTCCCCAGGT,C,C,C,Allows for saddle pattern to supercede tan points\n'
        '1,CoatPattern,agouti,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,G,G,G,Allows for saddle pattern to supercede tan point pattern\n'
        '2,CoatPattern,agouti,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,G,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,Allows for saddle pattern to supercede tan point pattern\n'
        '3,CoatPattern,agouti,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,G,G,G,Allows for saddle pattern to supercede tan point pattern\n'
        '4,CoatPattern,agouti,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,G,G,GTCCCCAGGTCAGAGTT,Allows for saddle pattern to supercede tan point pattern\n'
        '5,CoatPattern,agouti,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,G,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,Allows for saddle pattern to supercede tan point pattern\n'
        '1,CoatPattern,mask,Facial mask - variant p.(M264V),MC1R,T,C,C,T,Determines if a dog has a dark facial mask\n'
        '2,CoatPattern,mask,Facial mask - variant p.(M264V),MC1R,T,C,T,T,Determines if a dog has a dark facial mask\n'
        '3,CoatPattern,mask,Facial mask - variant p.(M264V),MC1R,T,C,C,C,Determines if a dog has a dark facial mask\n'
        '4,CoatPattern,mask,Facial mask - variant p.(M264V),MC1R,T,C,T,T,Determines if a dog has a dark facial mask\n'
        '5,CoatPattern,mask,Facial mask - variant p.(M264V),MC1R,T,C,C,T,Determines if a dog has a dark facial mask\n'
        '1,CoatPattern,grizzle,Sighthound grizzle - variant p.(G78V),MC1R,C,A,C,C,Determines if dog has shaded sable with widows peak pattern\n'
        '2,CoatPattern,grizzle,Sighthound grizzle - variant p.(G78V),MC1R,C,A,C,C,Determines if dog has shaded sable with widows peak pattern\n'
        '3,CoatPattern,grizzle,Sighthound grizzle - variant p.(G78V),MC1R,C,A,C,C,Determines if dog has shaded sable with widows peak pattern\n'
        '4,CoatPattern,grizzle,Sighthound grizzle - variant p.(G78V),MC1R,C,A,C,C,Determines if dog has shaded sable with widows peak pattern\n'
        '5,CoatPattern,grizzle,Sighthound grizzle - variant p.(G78V),MC1R,C,A,C,C,Determines if dog has shaded sable with widows peak pattern\n'
        '1,CoatPattern,northern_domino,Northern domino - variant p.(R301C),MC1R,G,A,G,G,Determines if a dog has domino facial and body pattern common in Northern breeds\n'
        '2,CoatPattern,northern_domino,Northern domino - variant p.(R301C),MC1R,G,A,G,G,Determines if a dog has domino facial and body pattern common in Northern breeds\n'
        '3,CoatPattern,northern_domino,Northern domino - variant p.(R301C),MC1R,G,A,G,G,Determines if a dog has domino facial and body pattern common in Northern breeds\n'
        '4,CoatPattern,northern_domino,Northern domino - variant p.(R301C),MC1R,G,A,G,G,Determines if a dog has domino facial and body pattern common in Northern breeds\n'
        '5,CoatPattern,northern_domino,Northern domino - variant p.(R301C),MC1R,G,A,G,G,Determines if a dog has domino facial and body pattern common in Northern breeds\n'
        '1,CoatPattern,red,Recessive red - variant p.(R306*),MC1R,G,A,G,A,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '2,CoatPattern,red,Recessive red - variant p.(R306*),MC1R,G,A,A,A,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '3,CoatPattern,red,Recessive red - variant p.(R306*),MC1R,G,A,G,G,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '4,CoatPattern,red,Recessive red - variant p.(R306*),MC1R,G,A,G,A,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '5,CoatPattern,red,Recessive red - variant p.(R306*),MC1R,G,A,G,G,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '1,CoatPattern,red,Recessive red - regulatory variant,MC1R,C,G,C,C,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '2,CoatPattern,red,Recessive red - regulatory variant,MC1R,C,G,C,C,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '3,CoatPattern,red,Recessive red - regulatory variant,MC1R,C,G,C,C,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '4,CoatPattern,red,Recessive red - regulatory variant,MC1R,C,G,C,C,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '5,CoatPattern,red,Recessive red - regulatory variant,MC1R,C,G,C,C,Pheomelanin solidly overrides eumelanin across entire coat\n'
        '1,CoatPattern,brindle,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC,TCCC,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '2,CoatPattern,brindle,Dominant black - variant p.(G23del),CBD103,TCCC,T,T,T,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '3,CoatPattern,brindle,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC,T,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '4,CoatPattern,brindle,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC,T,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '5,CoatPattern,brindle,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC,TCCC,Eumelanin solidly overrides pheomelanin across entire coat\n'
        '1,CoatPattern,brindle,Brindle - marker 1,intergenic,A,AGG,A,A,Associated with brindle coat pattern\n'
        '2,CoatPattern,brindle,Brindle - marker 1,intergenic,A,AGG,A,A,Associated with brindle coat pattern\n'
        '3,CoatPattern,brindle,Brindle - marker 1,intergenic,A,AGG,A,A,Associated with brindle coat pattern\n'
        '4,CoatPattern,brindle,Brindle - marker 1,intergenic,A,AGG,A,A,Associated with brindle coat pattern\n'
        '5,CoatPattern,brindle,Brindle - marker 1,intergenic,A,AGG,A,A,Associated with brindle coat pattern\n'
        '1,CoatPattern,brindle,Brindle - marker 2,intergenic,GCTTCCCTAAAA,G,GCTTCCCTAAAA,GCTTCCCTAAAA,Associated with brindle coat pattern\n'
        '2,CoatPattern,brindle,Brindle - marker 2,intergenic,GCTTCCCTAAAA,G,GCTTCCCTAAAA,GCTTCCCTAAAA,Associated with brindle coat pattern\n'
        '3,CoatPattern,brindle,Brindle - marker 2,intergenic,GCTTCCCTAAAA,G,GCTTCCCTAAAA,GCTTCCCTAAAA,Associated with brindle coat pattern\n'
        '4,CoatPattern,brindle,Brindle - marker 2,intergenic,GCTTCCCTAAAA,G,GCTTCCCTAAAA,GCTTCCCTAAAA,Associated with brindle coat pattern\n'
        '5,CoatPattern,brindle,Brindle - marker 2,intergenic,GCTTCCCTAAAA,G,GCTTCCCTAAAA,GCTTCCCTAAAA,Associated with brindle coat pattern\n'
        '1,CoatPattern,roan,Ticking - marker,USH2A,G,A,G,G,NA\n'
        '2,CoatPattern,roan,Ticking - marker,USH2A,G,A,G,G,NA\n'
        '3,CoatPattern,roan,Ticking - marker,USH2A,G,A,G,A,NA\n'
        '4,CoatPattern,roan,Ticking - marker,USH2A,G,A,G,G,NA\n'
        '5,CoatPattern,roan,Ticking - marker,USH2A,G,A,G,A,NA\n'
        '1,CoatPattern,merle,Harlequin - variant p.(V49I),PSMB7,T,G,T,T,Increases the amount of white on a merle patterned dog\n'
        '2,CoatPattern,merle,Harlequin - variant p.(V49I),PSMB7,T,G,T,T,Increases the amount of white on a merle patterned dog\n'
        '3,CoatPattern,merle,Harlequin - variant p.(V49I),PSMB7,T,G,T,T,Increases the amount of white on a merle patterned dog\n'
        '4,CoatPattern,merle,Harlequin - variant p.(V49I),PSMB7,T,G,T,T,Increases the amount of white on a merle patterned dog\n'
        '5,CoatPattern,merle,Harlequin - variant p.(V49I),PSMB7,T,G,T,T,Increases the amount of white on a merle patterned dog\n'
        '1,CoatType,coat_texture,Curly coat - variant p.(R151W),KRT71,C,T,C,C,Makes coat curly or wavy\n'
        '2,CoatType,coat_texture,Curly coat - variant p.(R151W),KRT71,C,T,C,C,Makes coat curly or wavy\n'
        '3,CoatType,coat_texture,Curly coat - variant p.(R151W),KRT71,C,T,C,T,Makes coat curly or wavy\n'
        '4,CoatType,coat_texture,Curly coat - variant p.(R151W),KRT71,C,T,C,C,Makes coat curly or wavy\n'
        '5,CoatType,coat_texture,Curly coat - variant p.(R151W),KRT71,C,T,C,C,Makes coat curly or wavy\n'
        '1,CoatType,coat_texture,Curly coat - variant p.(S422Rfs),KRT71,CTG,C,CTG,CTG,Makes coat curly or wavy\n'
        '2,CoatType,coat_texture,Curly coat - variant p.(S422Rfs),KRT71,CTG,C,CTG,CTG,Makes coat curly or wavy\n'
        '3,CoatType,coat_texture,Curly coat - variant p.(S422Rfs),KRT71,CTG,C,CTG,CTG,Makes coat curly or wavy\n'
        '4,CoatType,coat_texture,Curly coat - variant p.(S422Rfs),KRT71,CTG,C,CTG,CTG,Makes coat curly or wavy\n'
        '5,CoatType,coat_texture,Curly coat - variant p.(S422Rfs),KRT71,CTG,C,CTG,CTG,Makes coat curly or wavy\n'
        '1,CoatType,coat_length,Long coat - variant p.(C95F),FGF5,G,T,G,G,Increases coat length\n'
        '2,CoatType,coat_length,Long coat - variant p.(C95F),FGF5,G,T,G,G,Increases coat length\n'
        '3,CoatType,coat_length,Long coat - variant p.(C95F),FGF5,G,T,T,T,Increases coat length\n'
        '4,CoatType,coat_length,Long coat - variant p.(C95F),FGF5,G,T,G,G,Increases coat length\n'
        '5,CoatType,coat_length,Long coat - variant p.(C95F),FGF5,G,T,G,G,Increases coat length\n'
        '1,CoatType,furnishings,Furnishings - marker,RSPO2,A,C,A,A,Causes long hair on eyebrows and jowls\n'
        '2,CoatType,furnishings,Furnishings - marker,RSPO2,A,C,A,A,Causes long hair on eyebrows and jowls\n'
        '3,CoatType,furnishings,Furnishings - marker,RSPO2,A,C,A,A,Causes long hair on eyebrows and jowls\n'
        '4,CoatType,furnishings,Furnishings - marker,RSPO2,A,C,A,A,Causes long hair on eyebrows and jowls\n'
        '5,CoatType,furnishings,Furnishings - marker,RSPO2,A,C,A,A,Causes long hair on eyebrows and jowls\n'
        '1,CoatType,shedding,Shedding propensity - variant p.(A237T),MC5R,T,C,T,T,Decreases shedding\n'
        '2,CoatType,shedding,Shedding propensity - variant p.(A237T),MC5R,T,C,T,C,Decreases shedding\n'
        '3,CoatType,shedding,Shedding propensity - variant p.(A237T),MC5R,T,C,C,C,Decreases shedding\n'
        '4,CoatType,shedding,Shedding propensity - variant p.(A237T),MC5R,T,C,T,T,Decreases shedding\n'
        '5,CoatType,shedding,Shedding propensity - variant p.(A237T),MC5R,T,C,T,C,Decreases shedding\n'
        '1,CoatType,coat_layer,Single-layer coat - marker 1,ADRB1-AU1,C,T,T,T,Loss of double-layered coat\n'
        '2,CoatType,coat_layer,Single-layer coat - marker 1,ADRB1-AU1,C,T,T,T,Loss of double-layered coat\n'
        '3,CoatType,coat_layer,Single-layer coat - marker 1,ADRB1-AU1,C,T,C,C,Loss of double-layered coat\n'
        '4,CoatType,coat_layer,Single-layer coat - marker 1,ADRB1-AU1,C,T,T,T,Loss of double-layered coat\n'
        '5,CoatType,coat_layer,Single-layer coat - marker 1,ADRB1-AU1,C,T,T,C,Loss of double-layered coat\n'
        '1,CoatType,coat_layer,Single-layer coat - marker 2,ADRB1-AU1,G,A,A,A,Loss of double-layered coat\n'
        '2,CoatType,coat_layer,Single-layer coat - marker 2,ADRB1-AU1,G,A,G,G,Loss of double-layered coat\n'
        '3,CoatType,coat_layer,Single-layer coat - marker 2,ADRB1-AU1,G,A,G,G,Loss of double-layered coat\n'
        '4,CoatType,coat_layer,Single-layer coat - marker 2,ADRB1-AU1,G,A,A,A,Loss of double-layered coat\n'
        '5,CoatType,coat_layer,Single-layer coat - marker 2,ADRB1-AU1,G,A,A,G,Loss of double-layered coat\n'
        '1,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 1,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '2,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 1,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '3,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 1,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '4,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 1,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '5,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 1,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '1,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 2,EPAS1,G,T,T,T,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '2,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 2,EPAS1,G,T,T,T,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '3,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 2,EPAS1,G,T,T,T,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '4,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 2,EPAS1,G,T,T,T,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '5,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 2,EPAS1,G,T,T,T,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '1,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 3,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '2,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 3,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '3,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 3,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '4,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 3,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '5,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 3,EPAS1,G,A,G,G,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '1,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 4,EPAS1,C,T,C,C,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '2,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 4,EPAS1,C,T,C,C,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '3,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 4,EPAS1,C,T,C,C,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '4,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 4,EPAS1,C,T,C,C,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '5,SpecialFeatures,altitude,High altitude hypoxia tolerance - marker 4,EPAS1,C,T,C,C,"Observed in breeds from high elevations, like Tibetan mastiffs"\n'
        '1,SpecialFeatures,eye_color,Blue eyes - marker,ALX4,C,T,C,C,Associated with blue-colored eyes\n'
        '2,SpecialFeatures,eye_color,Blue eyes - marker,ALX4,C,T,C,T,Associated with blue-colored eyes\n'
        '3,SpecialFeatures,eye_color,Blue eyes - marker,ALX4,C,T,C,C,Associated with blue-colored eyes\n'
        '4,SpecialFeatures,eye_color,Blue eyes - marker,ALX4,C,T,C,C,Associated with blue-colored eyes\n'
        '5,SpecialFeatures,eye_color,Blue eyes - marker,ALX4,C,T,C,C,Associated with blue-colored eyes\n'
        '1,SpecialFeatures,leg_length,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,G,A,A,"Observed in breeds with short legs, such as Corgis"\n'
        '2,SpecialFeatures,leg_length,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,G,A,A,"Observed in breeds with short legs, such as Corgis"\n'
        '3,SpecialFeatures,leg_length,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,G,A,A,"Observed in breeds with short legs, such as Corgis"\n'
        '4,SpecialFeatures,leg_length,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,G,G,G,"Observed in breeds with short legs, such as Corgis"\n'
        '5,SpecialFeatures,leg_length,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,G,A,A,"Observed in breeds with short legs, such as Corgis"\n'
        '1,SpecialFeatures,leg_length,Long legs - marker 1,ESR1,C,T,T,T,Associated with longer than average legs\n'
        '2,SpecialFeatures,leg_length,Long legs - marker 1,ESR1,C,T,T,T,Associated with longer than average legs\n'
        '3,SpecialFeatures,leg_length,Long legs - marker 1,ESR1,C,T,T,T,Associated with longer than average legs\n'
        '4,SpecialFeatures,leg_length,Long legs - marker 1,ESR1,C,T,T,T,Associated with longer than average legs\n'
        '5,SpecialFeatures,leg_length,Long legs - marker 1,ESR1,C,T,T,T,Associated with longer than average legs\n'
        '1,SpecialFeatures,leg_length,Long legs - marker 2,ESR1,A,G,G,G,Associated with longer than average legs\n'
        '2,SpecialFeatures,leg_length,Long legs - marker 2,ESR1,A,G,G,G,Associated with longer than average legs\n'
        '3,SpecialFeatures,leg_length,Long legs - marker 2,ESR1,A,G,G,G,Associated with longer than average legs\n'
        '4,SpecialFeatures,leg_length,Long legs - marker 2,ESR1,A,G,G,G,Associated with longer than average legs\n'
        '5,SpecialFeatures,leg_length,Long legs - marker 2,ESR1,A,G,G,G,Associated with longer than average legs\n'
        '1,SpecialFeatures,leg_length,Long legs - marker 3,ESR1,T,G,G,G,Associated with longer than average legs\n'
        '2,SpecialFeatures,leg_length,Long legs - marker 3,ESR1,T,G,T,G,Associated with longer than average legs\n'
        '3,SpecialFeatures,leg_length,Long legs - marker 3,ESR1,T,G,G,G,Associated with longer than average legs\n'
        '4,SpecialFeatures,leg_length,Long legs - marker 3,ESR1,T,G,G,G,Associated with longer than average legs\n'
        '5,SpecialFeatures,leg_length,Long legs - marker 3,ESR1,T,G,G,G,Associated with longer than average legs\n'
        '1,SpecialFeatures,tail_length,Natural bob tail - variant p.(I63M),T,G,C,NA,NA,Dogs with variant born with shorter tail\n'
        '2,SpecialFeatures,tail_length,Natural bob tail - variant p.(I63M),T,G,C,G,G,Dogs with variant born with shorter tail\n'
        '3,SpecialFeatures,tail_length,Natural bob tail - variant p.(I63M),T,G,C,G,G,Dogs with variant born with shorter tail\n'
        '4,SpecialFeatures,tail_length,Natural bob tail - variant p.(I63M),T,G,C,G,G,Dogs with variant born with shorter tail\n'
        '5,SpecialFeatures,tail_length,Natural bob tail - variant p.(I63M),T,G,C,G,G,Dogs with variant born with shorter tail'
    )

def trail_pheno():
    return StringIO(
        'sample\ttab\tresult\tstring\n'
        '1\tCoatColor\tblack\tBlack\n'
        '1\tCoatColor\ttan\tTan\n'
        '1\tCoatPattern\tsable\tSable\n'
        '1\tCoatPattern\tmask\tFacial Mask\n'
        '1\tCoatPattern\tno_ticking\t\n'
        '1\tCoatPattern\tnot_brindle\t\n'
        '1\tCoatPattern\tnot_merle\t\n'
        '1\tCoatType\tstraight_coat\tStraight coat\n'
        '1\tCoatType\tshort_coat\tShort coat\n'
        '1\tCoatType\tno_furnishings\tNo eyebrow and muzzle furnishings\n'
        '1\tCoatType\tnormal_shedding\tNormal shedding\n'
        '1\tSpecialFeatures\tnormal_tail\tNormal length tail\n'
        '1\tSpecialFeatures\tnormal_legs\tNormal leg length\n'
        '1\tSpecialFeatures\tnot_hypoxia_adapted\tNo adaptation to high altitudes\n'
        '2\tCoatColor\tblack_leathers\tBlack (nose and paws)\n'
        '2\tCoatColor\tred\tRed\n'
        '2\tCoatPattern\ttan_points_hidden\tTan Points (hidden)\n'
        '2\tCoatPattern\trec_red\tNo mask, grizzle, or domino patterns\n'
        '2\tCoatPattern\tno_ticking\t\n'
        '2\tCoatPattern\tnot_brindle\t\n'
        '2\tCoatPattern\tnot_merle\t\n'
        '2\tCoatType\tstraight_coat\tStraight coat\n'
        '2\tCoatType\tshort_coat\tShort coat\n'
        '2\tCoatType\tno_furnishings\tNo eyebrow and muzzle furnishings\n'
        '2\tCoatType\tlow_shedding\tLow shedding\n'
        '2\tSpecialFeatures\tnormal_tail\tNormal length tail\n'
        '2\tSpecialFeatures\tnormal_legs\tNormal leg length\n'
        '2\tSpecialFeatures\tnot_hypoxia_adapted\tNo adaptation to high altitudes\n'
        '3\tCoatColor\tblack_solid\tBlack (solid coat)\n'
        '3\tCoatColor\ttan\tTan\n'
        '3\tCoatPattern\ttan_points_hidden\tTan Points (hidden)\n'
        '3\tCoatPattern\tmask_hidden\tFacial Mask (hidden)\n'
        '3\tCoatPattern\tno_ticking\t\n'
        '3\tCoatPattern\tnot_brindle\t\n'
        '3\tCoatPattern\tnot_merle\t\n'
        '3\tCoatType\twavy_coat\tWavy coat\n'
        '3\tCoatType\tlong_coat\tLong coat\n'
        '3\tCoatType\tno_furnishings\tNo eyebrow and muzzle furnishings\n'
        '3\tCoatType\tlow_shedding\tLow shedding\n'
        '3\tSpecialFeatures\tnormal_tail\tNormal length tail\n'
        '3\tSpecialFeatures\tnormal_legs\tNormal leg length\n'
        '3\tSpecialFeatures\tnot_hypoxia_adapted\tNo adaptation to high altitudes\n'
        '4\tCoatColor\tblack_solid\tBlack (solid coat)\n'
        '4\tCoatColor\tred\tRed\n'
        '4\tCoatPattern\tsable_hidden\tSable (hidden)\n'
        '4\tCoatPattern\tnormal_extension_hidden\tNo mask, grizzle, or domino patterns\n'
        '4\tCoatPattern\tno_ticking\t\n'
        '4\tCoatPattern\tnot_brindle\t\n'
        '4\tCoatPattern\tnot_merle\t\n'
        '4\tCoatType\tstraight_coat\tStraight coat\n'
        '4\tCoatType\tshort_coat\tShort coat\n'
        '4\tCoatType\tno_furnishings\tNo eyebrow and muzzle furnishings\n'
        '4\tCoatType\tnormal_shedding\tNormal shedding\n'
        '4\tSpecialFeatures\tnormal_tail\tNormal length tail\n'
        '4\tSpecialFeatures\tshort_legs_marker\tShortened leg length\n'
        '4\tSpecialFeatures\tnot_hypoxia_adapted\tNo adaptation to high altitudes\n'
        '5\tCoatColor\tblack\tBlack\n'
        '5\tCoatColor\ttan\tTan\n'
        '5\tCoatPattern\ttan_points\tTan Points\n'
        '5\tCoatPattern\tmask\tFacial Mask\n'
        '5\tCoatPattern\tno_ticking\t\n'
        '5\tCoatPattern\tnot_brindle\t\n'
        '5\tCoatPattern\tnot_merle\t\n'
        '5\tCoatType\tstraight_coat\tStraight coat\n'
        '5\tCoatType\tshort_coat\tShort coat\n'
        '5\tCoatType\tno_furnishings\tNo eyebrow and muzzle furnishings\n'
        '5\tCoatType\tlow_shedding\tLow shedding\n'
        '5\tSpecialFeatures\tnormal_tail\tNormal length tail\n'
        '5\tSpecialFeatures\tnormal_legs\tNormal leg length\n'
        '5\tSpecialFeatures\tnot_hypoxia_adapted\tNo adaptation to high altitudes'
    )

@pytest.fixture()
def old_outputs():
    return {
        'genotypeTable.csv': genotype_table(),
        'phenotypeTable.csv': phenotype_table(),
        'jsonTable.csv': json_table(),
        'trailblazerGenotypeTable.csv': trail_geno(),
        'trailblazerPhenotypeTable.tsv': trail_pheno(),
    }

@pytest.fixture()
def ids():
    return StringIO(
        '1\n'
        '22\n'
        '333\n'
    )


@pytest.fixture()
def roh_coi_csv():
    return StringIO(
        "id,nSeg,kbTot,kbAvg,coi\n"
        "1,1838,396380,215.658,0.179865\n"
        "22,2165,519333,239.876,0.235657\n"
        "333,1087,176535,162.406,0.0801061\n"
        "334,1088,176535,162.406,0.0801061\n"
    )


@pytest.fixture()
def report_pheno():
    return StringIO(
        'sample,tab,trait,result,string,value,color,image\n'
        '1,CoatColor,black,black,Black,,#17181D,\n'
        '22,CoatColor,black,black_leathers,Black (nose and paws),,#17181D,\n'
        '333,CoatColor,black,black_solid,Black (solid coat),,#17181D,\n'
        '1,CoatColor,red,tan,Tan,,#A86B39,\n'
        '22,CoatColor,red,red,Red,,#7E341B,\n'
        '333,CoatColor,red,tan,Tan,,#A86B39,\n'
        '1,CoatPattern,agouti,sable,Sable,,,CoatPattern_locusA_sable.svg\n'
        '22,CoatPattern,agouti,tan_points_hidden,Tan Points (hidden),,,CoatPattern_locusA_tanPoints.svg\n'
        '333,CoatPattern,agouti,tan_points_hidden,Tan Points (hidden),,,CoatPattern_locusA_tanPoints.svg\n'
        '1,CoatPattern,extension,mask,Facial Mask,,,CoatPattern_locusE_mask.svg\n'
        '22,CoatPattern,extension,rec_red,"No mask, grizzle, or domino patterns",,,\n'
        '333,CoatPattern,extension,mask_hidden,Facial Mask (hidden),,,CoatPattern_locusE_mask.svg\n'
        '1,CoatPattern,ticking,no_ticking,,,,\n'
        '22,CoatPattern,ticking,no_ticking,,,,\n'
        '333,CoatPattern,ticking,no_ticking,,,,\n'
        '1,CoatPattern,brindle,not_brindle,,,,\n'
        '22,CoatPattern,brindle,not_brindle,,,,\n'
        '333,CoatPattern,brindle,not_brindle,,,,\n'
        '1,CoatPattern,merle,not_merle,,,,\n'
        '22,CoatPattern,merle,not_merle,,,,\n'
        '333,CoatPattern,merle,not_merle,,,,\n'
        '1,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '22,CoatType,Coat Texture,straight_coat,Straight coat,,,CoatType_Curl_straight.svg\n'
        '333,CoatType,Coat Texture,wavy_coat,Wavy coat,,,CoatType_Curl_wavy.svg\n'
        '1,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '22,CoatType,Coat Length,short_coat,Short coat,,,CoatType_Length_short.svg\n'
        '333,CoatType,Coat Length,long_coat,Long coat,,,CoatType_Length_long.svg\n'
        '1,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '22,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '333,CoatType,Coat Furnishings,no_furnishings,No eyebrow and muzzle furnishings,,,CoatType_Furnishings_unfurnished.svg\n'
        '1,CoatType,Shedding Propensity,normal_shedding,Normal shedding,,,CoatType_Shedding_normal.svg\n'
        '22,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '333,CoatType,Shedding Propensity,low_shedding,Low shedding,,,CoatType_Shedding_low.svg\n'
        '1,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '22,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '333,SpecialFeatures,Skeletal - Tail Length,normal_tail,Normal length tail,,,SpecialFeatures_Tail_normal_tail.svg\n'
        '1,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '22,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '333,SpecialFeatures,Skeletal - Leg Length,normal_legs,Normal leg length,,,SpecialFeatures_Limbs_normal.svg\n'
        '1,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '22,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '333,SpecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
        '334,SPecialFeatures,High Altitude Adaptation,not_hypoxia_adapted,No adaptation to high altitudes,,,SpecialFeatures_Altitudes_not_adapted.svg\n'
    )


@pytest.fixture()
def report_json_table():
    return StringIO(
        'sample,tab,name,gene,firstCopy,secondCopy,possibleAlleles,effect\n'
        '1,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '333,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '22,CoatColor,Liver - variant p.(C41S),TYRP1,T,T,T & A,0\n'
        '1,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '333,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '22,CoatColor,Liver - variant p.(P345del),TYRP1,ACCT,ACCT,ACCT & A,0\n'
        '1,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '333,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '22,CoatColor,Liver - variant p.(Gln331*),TYRP1,C,C,C & T,0\n'
        '1,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '333,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '22,CoatColor,Cocoa - variant p.(T807*),HPS3,G,G,G & A,0\n'
        '1,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '333,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '22,CoatColor,Dilution - splice variant,MLPH,G,G,G & A,0\n'
        '1,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '333,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '22,CoatColor,Dilution - variant p.(Q235H),MLPH,G,G,G & C,0\n'
        '1,CoatColor,Red intensity - marker 1,lincRNA,T,A,T & A,0\n'
        '333,CoatColor,Red intensity - marker 1,lincRNA,T,A,T & A,0\n'
        '22,CoatColor,Red intensity - marker 1,lincRNA,A,A,T & A,0\n'
        '1,CoatColor,Red intensity - marker 2,intergenic,C,C,T & C,0\n'
        '333,CoatColor,Red intensity - marker 2,intergenic,T,C,T & C,0\n'
        '22,CoatColor,Red intensity - marker 2,intergenic,T,C,T & C,0\n'
        '1,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '333,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '22,CoatColor,Red intensity - marker 3,SLC264A,T,T,T & C,0\n'
        '1,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '333,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '22,CoatColor,Red intensity - marker 4,intergenic,C,C,T & C,0\n'
        '1,CoatColor,Red intensity - marker 5,TYR,G,G,G & A,0\n'
        '333,CoatColor,Red intensity - marker 5,TYR,G,A,G & A,0\n'
        '22,CoatColor,Red intensity - marker 5,TYR,G,A,G & A,0\n'
        '1,CoatPattern,Sable - variant p.(A82S),ASIP,T,T,G & T,0\n'
        '333,CoatPattern,Sable - variant p.(A82S),ASIP,G,G,G & T,0\n'
        '22,CoatPattern,Sable - variant p.(A82S),ASIP,G,G,G & T,0\n'
        '1,CoatPattern,Sable - variant p.(R83H),ASIP,A,A,G & A,0\n'
        '333,CoatPattern,Sable - variant p.(R83H),ASIP,G,G,G & A,0\n'
        '22,CoatPattern,Sable - variant p.(R83H),ASIP,G,G,G & A,0\n'
        '1,CoatPattern,Tan points - marker,ASIP,C,C,C & T,0\n'
        '333,CoatPattern,Tan points - marker,ASIP,T,T,C & T,0\n'
        '22,CoatPattern,Tan points - marker,ASIP,T,T,C & T,0\n'
        '1,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,C,C & T,0\n'
        '333,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,T,C & T,0\n'
        '22,CoatPattern,Recessive black - variant p.(R96C),ASIP,C,C,C & T,0\n'
        '1,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '333,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '22,CoatPattern,Saddle - marker 1,RALY,C,C,CAGAGTTTCCCCAGGT & C,0\n'
        '1,CoatPattern,Saddle - marker 2,RALY,G,G,GTCCCCAGGTCAGAGTT & G,0\n'
        '333,CoatPattern,Saddle - marker 2,RALY,G,G,GTCCCCAGGTCAGAGTT & G,0\n'
        '22,CoatPattern,Saddle - marker 2,RALY,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT,GTCCCCAGGTCAGAGTT & G,0\n'
        '1,CoatPattern,Facial mask - variant p.(M264V),MC1R,C,T,T & C,0\n'
        '333,CoatPattern,Facial mask - variant p.(M264V),MC1R,C,C,T & C,0\n'
        '22,CoatPattern,Facial mask - variant p.(M264V),MC1R,T,T,T & C,0\n'
        '1,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '333,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '22,CoatPattern,Sighthound grizzle - variant p.(G78V),MC1R,C,C,C & A,0\n'
        '1,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '333,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '22,CoatPattern,Northern domino - variant p.(R301C),MC1R,G,G,G & A,0\n'
        '1,CoatPattern,Recessive red - variant p.(R306*),MC1R,G,A,G & A,0\n'
        '333,CoatPattern,Recessive red - variant p.(R306*),MC1R,G,G,G & A,0\n'
        '22,CoatPattern,Recessive red - variant p.(R306*),MC1R,A,A,G & A,0\n'
        '1,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '333,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '22,CoatPattern,Recessive red - regulatory variant,MC1R,C,C,C & G,0\n'
        '1,CoatPattern,Dominant black - variant p.(G23del),CBD103,TCCC,TCCC,TCCC & T,0\n'
        '333,CoatPattern,Dominant black - variant p.(G23del),CBD103,TCCC,T,TCCC & T,0\n'
        '22,CoatPattern,Dominant black - variant p.(G23del),CBD103,T,T,TCCC & T,0\n'
        '1,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '333,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '22,CoatPattern,Brindle - marker 1,intergenic,A,A,A & AGG,0\n'
        '1,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '333,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '22,CoatPattern,Brindle - marker 2,intergenic,GCTTCCCTAAAA,GCTTCCCTAAAA,GCTTCCCTAAAA & G,0\n'
        '1,CoatPattern,Ticking - marker,USH2A,G,G,G & A,0\n'
        '333,CoatPattern,Ticking - marker,USH2A,G,A,G & A,0\n'
        '22,CoatPattern,Ticking - marker,USH2A,G,G,G & A,0\n'
        '1,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '333,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '22,CoatPattern,Harlequin - variant p.(V49I),PSMB7,T,T,T & G,0\n'
        '1,CoatType,Curly coat - variant p.(R151W),KRT71,C,C,C & T,0\n'
        '333,CoatType,Curly coat - variant p.(R151W),KRT71,C,T,C & T,0\n'
        '22,CoatType,Curly coat - variant p.(R151W),KRT71,C,C,C & T,0\n'
        '1,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '333,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '22,CoatType,Curly coat - variant p.(S422Rfs),KRT71,CTG,CTG,CTG & C,0\n'
        '1,CoatType,Long coat - variant p.(C95F),FGF5,G,G,G & T,0\n'
        '333,CoatType,Long coat - variant p.(C95F),FGF5,T,T,G & T,0\n'
        '22,CoatType,Long coat - variant p.(C95F),FGF5,G,G,G & T,0\n'
        '1,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '333,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '22,CoatType,Furnishings - marker,RSPO2,A,A,A & C,0\n'
        '1,CoatType,Shedding propensity - variant p.(A237T),MC5R,T,T,T & C,0\n'
        '333,CoatType,Shedding propensity - variant p.(A237T),MC5R,C,C,T & C,0\n'
        '22,CoatType,Shedding propensity - variant p.(A237T),MC5R,T,C,T & C,0\n'
        '1,CoatType,Single-layer coat - marker 1,ADRB1-AU1,T,T,C & T,0\n'
        '333,CoatType,Single-layer coat - marker 1,ADRB1-AU1,C,C,C & T,0\n'
        '22,CoatType,Single-layer coat - marker 1,ADRB1-AU1,T,T,C & T,0\n'
        '1,CoatType,Single-layer coat - marker 2,ADRB1-AU1,A,A,G & A,0\n'
        '333,CoatType,Single-layer coat - marker 2,ADRB1-AU1,G,G,G & A,0\n'
        '22,CoatType,Single-layer coat - marker 2,ADRB1-AU1,G,G,G & A,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '333,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '22,SpecialFeatures,High altitude hypoxia tolerance - marker 1,EPAS1,G,G,G & A,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '333,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '22,SpecialFeatures,High altitude hypoxia tolerance - marker 2,EPAS1,T,T,G & T,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '333,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '22,SpecialFeatures,High altitude hypoxia tolerance - marker 3,EPAS1,G,G,G & A,0\n'
        '1,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '333,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '22,SpecialFeatures,High altitude hypoxia tolerance - marker 4,EPAS1,C,C,C & T,0\n'
        '1,SpecialFeatures,Blue eyes - marker,ALX4,C,C,C & T,0\n'
        '333,SpecialFeatures,Blue eyes - marker,ALX4,C,C,C & T,0\n'
        '22,SpecialFeatures,Blue eyes - marker,ALX4,C,T,C & T,0\n'
        '1,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '333,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '22,SpecialFeatures,Shortened legs - marker,FGF4 retrogene on chromosome 18,A,A,A & G,0\n'
        '1,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '333,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '22,SpecialFeatures,Long legs - marker 1,ESR1,T,T,C & T,0\n'
        '1,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '333,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '22,SpecialFeatures,Long legs - marker 2,ESR1,G,G,A & G,0\n'
        '1,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '333,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '22,SpecialFeatures,Long legs - marker 3,ESR1,G,G,T & G,0\n'
        '1,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
        '333,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
        '22,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
        '23,SpecialFeatures,Natural bob tail - variant p.(I63M),T,G,G,G & C,0\n'
    )

@pytest.fixture
def report_body_geno():
    return StringIO(
        '"id","name","gene","possibleAlleles","firstCopy","secondCopy","effect"\n'
        '"1","Chromosome #3 Position: 91085576","near LCORL","A & G","G","A",0\n'
        '"22","Chromosome #3 Position: 91085576","near LCORL","A & G","G","G",-1\n'
        '"333","Chromosome #3 Position: 91085576","near LCORL","A & G","G","A",0\n'
        '"1","Chromosome #4 Position: 67040898","GHR","C & T","C","C",-1\n'
        '"22","Chromosome #4 Position: 67040898","GHR","C & T","C","C",-1\n'
        '"333","Chromosome #4 Position: 67040898","GHR","C & T","C","T",0\n'
        '"1","Chromosome #6 Position: 22864474","HS3ST2","A & G","G","G",-1\n'
        '"22","Chromosome #6 Position: 22864474","HS3ST2","A & G","A","G",0\n'
        '"333","Chromosome #6 Position: 22864474","HS3ST2","A & G","A","A",1\n'
        '"1","Chromosome #10 Position: 8356059","HMGA2","G & T","G","G",1\n'
        '"22","Chromosome #10 Position: 8356059","HMGA2","G & T","G","G",1\n'
        '"333","Chromosome #10 Position: 8356059","HMGA2","G & T","G","G",1\n'
        '"1","Chromosome #12 Position: 33792879","near OGFRL1","G & A","G","G",1\n'
        '"22","Chromosome #12 Position: 33792879","near OGFRL1","G & A","G","G",1\n'
        '"333","Chromosome #12 Position: 33792879","near OGFRL1","G & A","A","A",-1\n'
        '"1","Chromosome #15 Position: 41219654","IGF1","T & C","T","T",1\n'
        '"22","Chromosome #15 Position: 41219654","IGF1","T & C","C","C",-1\n'
        '"333","Chromosome #15 Position: 41219654","IGF1","T & C","C","T",0\n'
        '"1","Chromosome #17 Position: 36295546","ANAPC1","C & T","C","C",1\n'
        '"22","Chromosome #17 Position: 36295546","ANAPC1","C & T","C","C",1\n'
        '"333","Chromosome #17 Position: 36295546","ANAPC1","C & T","C","C",1\n'
        '"1","Chromosome #18 Position: 20428564","FGF4 retrogene","G & GC","G","G",1\n'
        '"22","Chromosome #18 Position: 20428564","FGF4 retrogene","G & GC","G","G",1\n'
        '"333","Chromosome #18 Position: 20428564","FGF4 retrogene","G & GC","G","G",1\n'
        '"1","Chromosome #26 Position: 12761780","near MED13L","G & A","G","G",-1\n'
        '"22","Chromosome #26 Position: 12761780","near MED13L","G & A","G","G",-1\n'
        '"333","Chromosome #26 Position: 12761780","near MED13L","G & A","G","A",0\n'
        '"1","Chromosome #32 Position: 5421641","non-coding","A & T","A","A",-1\n'
        '"22","Chromosome #32 Position: 5421641","non-coding","A & T","A","A",-1\n'
        '"333","Chromosome #32 Position: 5421641","non-coding","A & T","A","A",-1\n'
        '"334","Chromosome #32 Position: 5421641","non-coding","A & T","A","A",-1\n'
    )

@pytest.fixture
def report_body_pheno():
    return StringIO(
        'id,prediction\n'
        '1,2.53063333333333\n'
        '22,2.38136666666667\n'
        '333,1.90556666666667\n'
        '334,1.90556666666663\n'
    )

@pytest.fixture
def report_white_geno():
    return StringIO(
        '"id","name","gene","possibleAlleles","firstCopy","secondCopy","effect"\n'
        '"1","Chromosome #4 Position: 4882111","non-coding","G & T","G","T",0\n'
        '"22","Chromosome #4 Position: 4882111","non-coding","G & T","G","G",1\n'
        '"333","Chromosome #4 Position: 4882111","non-coding","G & T","T","T",-1\n'
        '"1","Chromosome #5 Position: 63694334","MC1R","G & A","G","A",0\n'
        '"22","Chromosome #5 Position: 63694334","MC1R","G & A","A","A",-1\n'
        '"333","Chromosome #5 Position: 63694334","MC1R","G & A","G","G",1\n'
        '"1","Chromosome #14 Position: 29948181","AGMO","G & A","G","G",1\n'
        '"22","Chromosome #14 Position: 29948181","AGMO","G & A","G","G",1\n'
        '"333","Chromosome #14 Position: 29948181","AGMO","G & A","A","A",-1\n'
        '"1","Chromosome #20 Position: 21792546","MITF","G & A","G","G",1\n'
        '"22","Chromosome #20 Position: 21792546","MITF","G & A","G","G",1\n'
        '"333","Chromosome #20 Position: 21792546","MITF","G & A","G","G",1\n'
        '"1","Chromosome #20 Position: 21797796","MITF","A & C","A","A",1\n'
        '"22","Chromosome #20 Position: 21797796","MITF","A & C","A","A",1\n'
        '"333","Chromosome #20 Position: 21797796","MITF","A & C","C","C",-1\n'
        '"1","Chromosome #20 Position: 21825467","MITF","A & C","A","A",1\n'
        '"22","Chromosome #20 Position: 21825467","MITF","A & C","A","A",1\n'
        '"333","Chromosome #20 Position: 21825467","MITF","A & C","A","A",1\n'
        '"1","Chromosome #20 Position: 21827584","MITF","TTTTTTC & TTTTTTCTTTTTC","TTTTTTC","TTTTTTC",-1\n'
        '"22","Chromosome #20 Position: 21827584","MITF","TTTTTTC & TTTTTTCTTTTTC","TTTTTTC","TTTTTTC",-1\n'
        '"333","Chromosome #20 Position: 21827584","MITF","TTTTTTC & TTTTTTCTTTTTC","TTTTTTC","TTTTTTC",-1\n'
        '"1","Chromosome #20 Position: 21827657","MITF","T & TTTCTTTTC","T","T",1\n'
        '"22","Chromosome #20 Position: 21827657","MITF","T & TTTCTTTTC","T","T",1\n'
        '"333","Chromosome #20 Position: 21827657","MITF","T & TTTCTTTTC","TTTCTTTTC","TTTCTTTTC",-1\n'
        '"1","Chromosome #20 Position: 21829531","MITF","T & TA","T","T",1\n'
        '"22","Chromosome #20 Position: 21829531","MITF","T & TA","T","T",1\n'
        '"333","Chromosome #20 Position: 21829531","MITF","T & TA","TA","TA",-1\n'
        '"1","Chromosome #20 Position: 21834982","MITF","T & A","A","A",1\n'
        '"22","Chromosome #20 Position: 21834982","MITF","T & A","A","A",1\n'
        '"333","Chromosome #20 Position: 21834982","MITF","T & A","T","T",-1\n'
        '"334","Chromosome #20 Position: 21834982","MITF","T & A","T","T",-1\n'
    )

@pytest.fixture
def report_white_pheno():
    return StringIO(
        'id,prediction\n'
        '1,5\n'
        '22,5\n'
        '333,2\n'
        '334,3\n'
    )

@pytest.fixture
def get_adm():
    @contextmanager
    def adm_function(id):
        breed = id[0]
        yield StringIO(
            f'[{{"breed":{breed},"percent":0.9269}},{{"breed":245,"percent":0.0447}}]'
        )

    return adm_function

@pytest.fixture
def json_1():
    return StringIO(
        '''
{
    "id": "1",
    "inbreeding": 0.179865,
    "data": [
        {
            "breed": 1,
            "percent": 0.9269
        },
        {
            "breed": 245,
            "percent": 0.0447
        }
    ],
    "panels": [
        {
            "name": "Body Size",
            "id": "body-size",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    2.6064380952380914
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Chromosome #3 Position: 91085576",
                        "gene": "near LCORL",
                        "possibleAlleles": "A & G",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #4 Position: 67040898",
                        "gene": "GHR",
                        "possibleAlleles": "C & T",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #6 Position: 22864474",
                        "gene": "HS3ST2",
                        "possibleAlleles": "A & G",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #10 Position: 8356059",
                        "gene": "HMGA2",
                        "possibleAlleles": "G & T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #12 Position: 33792879",
                        "gene": "near OGFRL1",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #15 Position: 41219654",
                        "gene": "IGF1",
                        "possibleAlleles": "T & C",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #17 Position: 36295546",
                        "gene": "ANAPC1",
                        "possibleAlleles": "C & T",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #18 Position: 20428564",
                        "gene": "FGF4 retrogene",
                        "possibleAlleles": "G & GC",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #26 Position: 12761780",
                        "gene": "near MED13L",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #32 Position: 5421641",
                        "gene": "non-coding",
                        "possibleAlleles": "A & T",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": -1
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Colors",
            "id": "coat-color",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Black",
                    "#17181D",
                    "Tan",
                    "#A86B39"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Liver - variant p.(C41S)",
                        "gene": "TYRP1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & A",
                        "effect": 0
                    },
                    {
                        "name": "Liver - variant p.(P345del)",
                        "gene": "TYRP1",
                        "firstCopy": "ACCT",
                        "secondCopy": "ACCT",
                        "possibleAlleles": "ACCT & A",
                        "effect": 0
                    },
                    {
                        "name": "Liver - variant p.(Gln331*)",
                        "gene": "TYRP1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Cocoa - variant p.(T807*)",
                        "gene": "HPS3",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Dilution - splice variant",
                        "gene": "MLPH",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Dilution - variant p.(Q235H)",
                        "gene": "MLPH",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 1",
                        "gene": "lincRNA",
                        "firstCopy": "T",
                        "secondCopy": "A",
                        "possibleAlleles": "T & A",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 2",
                        "gene": "intergenic",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 3",
                        "gene": "SLC264A",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 4",
                        "gene": "intergenic",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 5",
                        "gene": "TYR",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Coat Pattern",
            "id": "coat-pattern",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "CoatPattern_locusA_sable.svg",
                    "Sable",
                    "CoatPattern_locusE_mask.svg",
                    "Facial Mask"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Sable - variant p.(A82S)",
                        "gene": "ASIP",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "Sable - variant p.(R83H)",
                        "gene": "ASIP",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Tan points - marker",
                        "gene": "ASIP",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Recessive black - variant p.(R96C)",
                        "gene": "ASIP",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Saddle - marker 1",
                        "gene": "RALY",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "CAGAGTTTCCCCAGGT & C",
                        "effect": 0
                    },
                    {
                        "name": "Saddle - marker 2",
                        "gene": "RALY",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "GTCCCCAGGTCAGAGTT & G",
                        "effect": 0
                    },
                    {
                        "name": "Facial mask - variant p.(M264V)",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "T",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Sighthound grizzle - variant p.(G78V)",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & A",
                        "effect": 0
                    },
                    {
                        "name": "Northern domino - variant p.(R301C)",
                        "gene": "MC1R",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Recessive red - variant p.(R306*)",
                        "gene": "MC1R",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Recessive red - regulatory variant",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & G",
                        "effect": 0
                    },
                    {
                        "name": "Dominant black - variant p.(G23del)",
                        "gene": "CBD103",
                        "firstCopy": "TCCC",
                        "secondCopy": "TCCC",
                        "possibleAlleles": "TCCC & T",
                        "effect": 0
                    },
                    {
                        "name": "Brindle - marker 1",
                        "gene": "intergenic",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & AGG",
                        "effect": 0
                    },
                    {
                        "name": "Brindle - marker 2",
                        "gene": "intergenic",
                        "firstCopy": "GCTTCCCTAAAA",
                        "secondCopy": "GCTTCCCTAAAA",
                        "possibleAlleles": "GCTTCCCTAAAA & G",
                        "effect": 0
                    },
                    {
                        "name": "Ticking - marker",
                        "gene": "USH2A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Harlequin - variant p.(V49I)",
                        "gene": "PSMB7",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & G",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "White Spotting",
            "id": "white-spotting",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    5
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Chromosome #4 Position: 4882111",
                        "gene": "non-coding",
                        "possibleAlleles": "G & T",
                        "firstCopy": "G",
                        "secondCopy": "T",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #5 Position: 63694334",
                        "gene": "MC1R",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #14 Position: 29948181",
                        "gene": "AGMO",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21792546",
                        "gene": "MITF",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21797796",
                        "gene": "MITF",
                        "possibleAlleles": "A & C",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21825467",
                        "gene": "MITF",
                        "possibleAlleles": "A & C",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21827584",
                        "gene": "MITF",
                        "possibleAlleles": "TTTTTTC & TTTTTTCTTTTTC",
                        "firstCopy": "TTTTTTC",
                        "secondCopy": "TTTTTTC",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21827657",
                        "gene": "MITF",
                        "possibleAlleles": "T & TTTCTTTTC",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21829531",
                        "gene": "MITF",
                        "possibleAlleles": "T & TA",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21834982",
                        "gene": "MITF",
                        "possibleAlleles": "T & A",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Coat Type",
            "id": "coat-type",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Coat Texture",
                    "CoatType_Curl_straight.svg",
                    "Straight coat",
                    "Coat Length",
                    "CoatType_Length_short.svg",
                    "Short coat",
                    "Coat Furnishings",
                    "CoatType_Furnishings_unfurnished.svg",
                    "No eyebrow and muzzle furnishings",
                    "Shedding Propensity",
                    "CoatType_Shedding_normal.svg",
                    "Normal shedding"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Curly coat - variant p.(R151W)",
                        "gene": "KRT71",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Curly coat - variant p.(S422Rfs)",
                        "gene": "KRT71",
                        "firstCopy": "CTG",
                        "secondCopy": "CTG",
                        "possibleAlleles": "CTG & C",
                        "effect": 0
                    },
                    {
                        "name": "Long coat - variant p.(C95F)",
                        "gene": "FGF5",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "Furnishings - marker",
                        "gene": "RSPO2",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & C",
                        "effect": 0
                    },
                    {
                        "name": "Shedding propensity - variant p.(A237T)",
                        "gene": "MC5R",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Single-layer coat - marker 1",
                        "gene": "ADRB1-AU1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Single-layer coat - marker 2",
                        "gene": "ADRB1-AU1",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Special Features",
            "id": "special-features",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Skeletal - Tail Length",
                    "SpecialFeatures_Tail_normal_tail.svg",
                    "Normal length tail",
                    "Skeletal - Leg Length",
                    "SpecialFeatures_Limbs_normal.svg",
                    "Normal leg length",
                    "High Altitude Adaptation",
                    "SpecialFeatures_Altitudes_not_adapted.svg",
                    "No adaptation to high altitudes"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "High altitude hypoxia tolerance - marker 1",
                        "gene": "EPAS1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 2",
                        "gene": "EPAS1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 3",
                        "gene": "EPAS1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 4",
                        "gene": "EPAS1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Blue eyes - marker",
                        "gene": "ALX4",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Shortened legs - marker",
                        "gene": "FGF4 retrogene on chromosome 18",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & G",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 1",
                        "gene": "ESR1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 2",
                        "gene": "ESR1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "A & G",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 3",
                        "gene": "ESR1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "T & G",
                        "effect": 0
                    },
                    {
                        "name": "Natural bob tail - variant p.(I63M)",
                        "gene": "T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & C",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        }
    ]
}
        '''
    )

@pytest.fixture
def json_22():
    return StringIO(
'''
{
    "id": "22",
    "inbreeding": 0.235657,
    "data": [
        {
            "breed": 2,
            "percent": 0.9269
        },
        {
            "breed": 245,
            "percent": 0.0447
        }
    ],
    "panels": [
        {
            "name": "Body Size",
            "id": "body-size",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    2.4358476190476233
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Chromosome #3 Position: 91085576",
                        "gene": "near LCORL",
                        "possibleAlleles": "A & G",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #4 Position: 67040898",
                        "gene": "GHR",
                        "possibleAlleles": "C & T",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #6 Position: 22864474",
                        "gene": "HS3ST2",
                        "possibleAlleles": "A & G",
                        "firstCopy": "A",
                        "secondCopy": "G",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #10 Position: 8356059",
                        "gene": "HMGA2",
                        "possibleAlleles": "G & T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #12 Position: 33792879",
                        "gene": "near OGFRL1",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #15 Position: 41219654",
                        "gene": "IGF1",
                        "possibleAlleles": "T & C",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #17 Position: 36295546",
                        "gene": "ANAPC1",
                        "possibleAlleles": "C & T",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #18 Position: 20428564",
                        "gene": "FGF4 retrogene",
                        "possibleAlleles": "G & GC",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #26 Position: 12761780",
                        "gene": "near MED13L",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #32 Position: 5421641",
                        "gene": "non-coding",
                        "possibleAlleles": "A & T",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": -1
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Colors",
            "id": "coat-color",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Black (nose and paws)",
                    "#17181D",
                    "Red",
                    "#7E341B"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Liver - variant p.(C41S)",
                        "gene": "TYRP1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & A",
                        "effect": 0
                    },
                    {
                        "name": "Liver - variant p.(P345del)",
                        "gene": "TYRP1",
                        "firstCopy": "ACCT",
                        "secondCopy": "ACCT",
                        "possibleAlleles": "ACCT & A",
                        "effect": 0
                    },
                    {
                        "name": "Liver - variant p.(Gln331*)",
                        "gene": "TYRP1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Cocoa - variant p.(T807*)",
                        "gene": "HPS3",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Dilution - splice variant",
                        "gene": "MLPH",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Dilution - variant p.(Q235H)",
                        "gene": "MLPH",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 1",
                        "gene": "lincRNA",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "T & A",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 2",
                        "gene": "intergenic",
                        "firstCopy": "T",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 3",
                        "gene": "SLC264A",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 4",
                        "gene": "intergenic",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 5",
                        "gene": "TYR",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Coat Pattern",
            "id": "coat-pattern",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "CoatPattern_locusA_tanPoints.svg",
                    "Tan Points (hidden)"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Sable - variant p.(A82S)",
                        "gene": "ASIP",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "Sable - variant p.(R83H)",
                        "gene": "ASIP",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Tan points - marker",
                        "gene": "ASIP",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Recessive black - variant p.(R96C)",
                        "gene": "ASIP",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Saddle - marker 1",
                        "gene": "RALY",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "CAGAGTTTCCCCAGGT & C",
                        "effect": 0
                    },
                    {
                        "name": "Saddle - marker 2",
                        "gene": "RALY",
                        "firstCopy": "GTCCCCAGGTCAGAGTT",
                        "secondCopy": "GTCCCCAGGTCAGAGTT",
                        "possibleAlleles": "GTCCCCAGGTCAGAGTT & G",
                        "effect": 0
                    },
                    {
                        "name": "Facial mask - variant p.(M264V)",
                        "gene": "MC1R",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Sighthound grizzle - variant p.(G78V)",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & A",
                        "effect": 0
                    },
                    {
                        "name": "Northern domino - variant p.(R301C)",
                        "gene": "MC1R",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Recessive red - variant p.(R306*)",
                        "gene": "MC1R",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Recessive red - regulatory variant",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & G",
                        "effect": 0
                    },
                    {
                        "name": "Dominant black - variant p.(G23del)",
                        "gene": "CBD103",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "TCCC & T",
                        "effect": 0
                    },
                    {
                        "name": "Brindle - marker 1",
                        "gene": "intergenic",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & AGG",
                        "effect": 0
                    },
                    {
                        "name": "Brindle - marker 2",
                        "gene": "intergenic",
                        "firstCopy": "GCTTCCCTAAAA",
                        "secondCopy": "GCTTCCCTAAAA",
                        "possibleAlleles": "GCTTCCCTAAAA & G",
                        "effect": 0
                    },
                    {
                        "name": "Ticking - marker",
                        "gene": "USH2A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Harlequin - variant p.(V49I)",
                        "gene": "PSMB7",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & G",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "White Spotting",
            "id": "white-spotting",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    5
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Chromosome #4 Position: 4882111",
                        "gene": "non-coding",
                        "possibleAlleles": "G & T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #5 Position: 63694334",
                        "gene": "MC1R",
                        "possibleAlleles": "G & A",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #14 Position: 29948181",
                        "gene": "AGMO",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21792546",
                        "gene": "MITF",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21797796",
                        "gene": "MITF",
                        "possibleAlleles": "A & C",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21825467",
                        "gene": "MITF",
                        "possibleAlleles": "A & C",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21827584",
                        "gene": "MITF",
                        "possibleAlleles": "TTTTTTC & TTTTTTCTTTTTC",
                        "firstCopy": "TTTTTTC",
                        "secondCopy": "TTTTTTC",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21827657",
                        "gene": "MITF",
                        "possibleAlleles": "T & TTTCTTTTC",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21829531",
                        "gene": "MITF",
                        "possibleAlleles": "T & TA",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21834982",
                        "gene": "MITF",
                        "possibleAlleles": "T & A",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Coat Type",
            "id": "coat-type",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Coat Texture",
                    "CoatType_Curl_straight.svg",
                    "Straight coat",
                    "Coat Length",
                    "CoatType_Length_short.svg",
                    "Short coat",
                    "Coat Furnishings",
                    "CoatType_Furnishings_unfurnished.svg",
                    "No eyebrow and muzzle furnishings",
                    "Shedding Propensity",
                    "CoatType_Shedding_low.svg",
                    "Low shedding"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Curly coat - variant p.(R151W)",
                        "gene": "KRT71",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Curly coat - variant p.(S422Rfs)",
                        "gene": "KRT71",
                        "firstCopy": "CTG",
                        "secondCopy": "CTG",
                        "possibleAlleles": "CTG & C",
                        "effect": 0
                    },
                    {
                        "name": "Long coat - variant p.(C95F)",
                        "gene": "FGF5",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "Furnishings - marker",
                        "gene": "RSPO2",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & C",
                        "effect": 0
                    },
                    {
                        "name": "Shedding propensity - variant p.(A237T)",
                        "gene": "MC5R",
                        "firstCopy": "T",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Single-layer coat - marker 1",
                        "gene": "ADRB1-AU1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Single-layer coat - marker 2",
                        "gene": "ADRB1-AU1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Special Features",
            "id": "special-features",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Skeletal - Tail Length",
                    "SpecialFeatures_Tail_normal_tail.svg",
                    "Normal length tail",
                    "Skeletal - Leg Length",
                    "SpecialFeatures_Limbs_normal.svg",
                    "Normal leg length",
                    "High Altitude Adaptation",
                    "SpecialFeatures_Altitudes_not_adapted.svg",
                    "No adaptation to high altitudes"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "High altitude hypoxia tolerance - marker 1",
                        "gene": "EPAS1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 2",
                        "gene": "EPAS1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 3",
                        "gene": "EPAS1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 4",
                        "gene": "EPAS1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Blue eyes - marker",
                        "gene": "ALX4",
                        "firstCopy": "C",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Shortened legs - marker",
                        "gene": "FGF4 retrogene on chromosome 18",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & G",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 1",
                        "gene": "ESR1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 2",
                        "gene": "ESR1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "A & G",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 3",
                        "gene": "ESR1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "T & G",
                        "effect": 0
                    },
                    {
                        "name": "Natural bob tail - variant p.(I63M)",
                        "gene": "T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & C",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        }
    ]
}
    '''
    )

@pytest.fixture
def json_333():
    return StringIO(
'''
{
    "id": "333",
    "inbreeding": 0.0801061,
    "data": [
        {
            "breed": 3,
            "percent": 0.9269
        },
        {
            "breed": 245,
            "percent": 0.0447
        }
    ],
    "panels": [
        {
            "name": "Body Size",
            "id": "body-size",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    1.8920761904761942
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Chromosome #3 Position: 91085576",
                        "gene": "near LCORL",
                        "possibleAlleles": "A & G",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #4 Position: 67040898",
                        "gene": "GHR",
                        "possibleAlleles": "C & T",
                        "firstCopy": "C",
                        "secondCopy": "T",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #6 Position: 22864474",
                        "gene": "HS3ST2",
                        "possibleAlleles": "A & G",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #10 Position: 8356059",
                        "gene": "HMGA2",
                        "possibleAlleles": "G & T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #12 Position: 33792879",
                        "gene": "near OGFRL1",
                        "possibleAlleles": "G & A",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #15 Position: 41219654",
                        "gene": "IGF1",
                        "possibleAlleles": "T & C",
                        "firstCopy": "C",
                        "secondCopy": "T",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #17 Position: 36295546",
                        "gene": "ANAPC1",
                        "possibleAlleles": "C & T",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #18 Position: 20428564",
                        "gene": "FGF4 retrogene",
                        "possibleAlleles": "G & GC",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #26 Position: 12761780",
                        "gene": "near MED13L",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "effect": 0
                    },
                    {
                        "name": "Chromosome #32 Position: 5421641",
                        "gene": "non-coding",
                        "possibleAlleles": "A & T",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": -1
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Colors",
            "id": "coat-color",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Black (solid coat)",
                    "#17181D",
                    "Tan",
                    "#A86B39"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Liver - variant p.(C41S)",
                        "gene": "TYRP1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & A",
                        "effect": 0
                    },
                    {
                        "name": "Liver - variant p.(P345del)",
                        "gene": "TYRP1",
                        "firstCopy": "ACCT",
                        "secondCopy": "ACCT",
                        "possibleAlleles": "ACCT & A",
                        "effect": 0
                    },
                    {
                        "name": "Liver - variant p.(Gln331*)",
                        "gene": "TYRP1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Cocoa - variant p.(T807*)",
                        "gene": "HPS3",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Dilution - splice variant",
                        "gene": "MLPH",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Dilution - variant p.(Q235H)",
                        "gene": "MLPH",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 1",
                        "gene": "lincRNA",
                        "firstCopy": "T",
                        "secondCopy": "A",
                        "possibleAlleles": "T & A",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 2",
                        "gene": "intergenic",
                        "firstCopy": "T",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 3",
                        "gene": "SLC264A",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 4",
                        "gene": "intergenic",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Red intensity - marker 5",
                        "gene": "TYR",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Coat Pattern",
            "id": "coat-pattern",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "CoatPattern_locusA_tanPoints.svg",
                    "Tan Points (hidden)",
                    "CoatPattern_locusE_mask.svg",
                    "Facial Mask (hidden)"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Sable - variant p.(A82S)",
                        "gene": "ASIP",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "Sable - variant p.(R83H)",
                        "gene": "ASIP",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Tan points - marker",
                        "gene": "ASIP",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Recessive black - variant p.(R96C)",
                        "gene": "ASIP",
                        "firstCopy": "C",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Saddle - marker 1",
                        "gene": "RALY",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "CAGAGTTTCCCCAGGT & C",
                        "effect": 0
                    },
                    {
                        "name": "Saddle - marker 2",
                        "gene": "RALY",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "GTCCCCAGGTCAGAGTT & G",
                        "effect": 0
                    },
                    {
                        "name": "Facial mask - variant p.(M264V)",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Sighthound grizzle - variant p.(G78V)",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & A",
                        "effect": 0
                    },
                    {
                        "name": "Northern domino - variant p.(R301C)",
                        "gene": "MC1R",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Recessive red - variant p.(R306*)",
                        "gene": "MC1R",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Recessive red - regulatory variant",
                        "gene": "MC1R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & G",
                        "effect": 0
                    },
                    {
                        "name": "Dominant black - variant p.(G23del)",
                        "gene": "CBD103",
                        "firstCopy": "TCCC",
                        "secondCopy": "T",
                        "possibleAlleles": "TCCC & T",
                        "effect": 0
                    },
                    {
                        "name": "Brindle - marker 1",
                        "gene": "intergenic",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & AGG",
                        "effect": 0
                    },
                    {
                        "name": "Brindle - marker 2",
                        "gene": "intergenic",
                        "firstCopy": "GCTTCCCTAAAA",
                        "secondCopy": "GCTTCCCTAAAA",
                        "possibleAlleles": "GCTTCCCTAAAA & G",
                        "effect": 0
                    },
                    {
                        "name": "Ticking - marker",
                        "gene": "USH2A",
                        "firstCopy": "G",
                        "secondCopy": "A",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "Harlequin - variant p.(V49I)",
                        "gene": "PSMB7",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "T & G",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "White Spotting",
            "id": "white-spotting",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    2
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Chromosome #4 Position: 4882111",
                        "gene": "non-coding",
                        "possibleAlleles": "G & T",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #5 Position: 63694334",
                        "gene": "MC1R",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #14 Position: 29948181",
                        "gene": "AGMO",
                        "possibleAlleles": "G & A",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21792546",
                        "gene": "MITF",
                        "possibleAlleles": "G & A",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21797796",
                        "gene": "MITF",
                        "possibleAlleles": "A & C",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21825467",
                        "gene": "MITF",
                        "possibleAlleles": "A & C",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "effect": 1
                    },
                    {
                        "name": "Chromosome #20 Position: 21827584",
                        "gene": "MITF",
                        "possibleAlleles": "TTTTTTC & TTTTTTCTTTTTC",
                        "firstCopy": "TTTTTTC",
                        "secondCopy": "TTTTTTC",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21827657",
                        "gene": "MITF",
                        "possibleAlleles": "T & TTTCTTTTC",
                        "firstCopy": "TTTCTTTTC",
                        "secondCopy": "TTTCTTTTC",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21829531",
                        "gene": "MITF",
                        "possibleAlleles": "T & TA",
                        "firstCopy": "TA",
                        "secondCopy": "TA",
                        "effect": -1
                    },
                    {
                        "name": "Chromosome #20 Position: 21834982",
                        "gene": "MITF",
                        "possibleAlleles": "T & A",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "effect": -1
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Coat Type",
            "id": "coat-type",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Coat Texture",
                    "CoatType_Curl_wavy.svg",
                    "Wavy coat",
                    "Coat Length",
                    "CoatType_Length_long.svg",
                    "Long coat",
                    "Coat Furnishings",
                    "CoatType_Furnishings_unfurnished.svg",
                    "No eyebrow and muzzle furnishings",
                    "Shedding Propensity",
                    "CoatType_Shedding_low.svg",
                    "Low shedding"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "Curly coat - variant p.(R151W)",
                        "gene": "KRT71",
                        "firstCopy": "C",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Curly coat - variant p.(S422Rfs)",
                        "gene": "KRT71",
                        "firstCopy": "CTG",
                        "secondCopy": "CTG",
                        "possibleAlleles": "CTG & C",
                        "effect": 0
                    },
                    {
                        "name": "Long coat - variant p.(C95F)",
                        "gene": "FGF5",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "Furnishings - marker",
                        "gene": "RSPO2",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & C",
                        "effect": 0
                    },
                    {
                        "name": "Shedding propensity - variant p.(A237T)",
                        "gene": "MC5R",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "T & C",
                        "effect": 0
                    },
                    {
                        "name": "Single-layer coat - marker 1",
                        "gene": "ADRB1-AU1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Single-layer coat - marker 2",
                        "gene": "ADRB1-AU1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        },
        {
            "name": "Special Features",
            "id": "special-features",
            "top": {
                "title": null,
                "content": "null",
                "genericValues": [
                    "Skeletal - Tail Length",
                    "SpecialFeatures_Tail_normal_tail.svg",
                    "Normal length tail",
                    "Skeletal - Leg Length",
                    "SpecialFeatures_Limbs_normal.svg",
                    "Normal leg length",
                    "High Altitude Adaptation",
                    "SpecialFeatures_Altitudes_not_adapted.svg",
                    "No adaptation to high altitudes"
                ]
            },
            "results": {
                "title": null,
                "content": "null",
                "traits": [
                    {
                        "name": "High altitude hypoxia tolerance - marker 1",
                        "gene": "EPAS1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 2",
                        "gene": "EPAS1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "G & T",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 3",
                        "gene": "EPAS1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & A",
                        "effect": 0
                    },
                    {
                        "name": "High altitude hypoxia tolerance - marker 4",
                        "gene": "EPAS1",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Blue eyes - marker",
                        "gene": "ALX4",
                        "firstCopy": "C",
                        "secondCopy": "C",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Shortened legs - marker",
                        "gene": "FGF4 retrogene on chromosome 18",
                        "firstCopy": "A",
                        "secondCopy": "A",
                        "possibleAlleles": "A & G",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 1",
                        "gene": "ESR1",
                        "firstCopy": "T",
                        "secondCopy": "T",
                        "possibleAlleles": "C & T",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 2",
                        "gene": "ESR1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "A & G",
                        "effect": 0
                    },
                    {
                        "name": "Long legs - marker 3",
                        "gene": "ESR1",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "T & G",
                        "effect": 0
                    },
                    {
                        "name": "Natural bob tail - variant p.(I63M)",
                        "gene": "T",
                        "firstCopy": "G",
                        "secondCopy": "G",
                        "possibleAlleles": "G & C",
                        "effect": 0
                    }
                ],
                "genericValues": null
            }
        }
    ]
}
    '''
    )
