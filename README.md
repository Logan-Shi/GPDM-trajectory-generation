This is a version of the GPDM code that you are free to use for academic purposes.  It is obviously not 'production code', and I do not guarantee its correctness or efficiency.  

See example.m for the code to learn a single walker model and generate samples from the learned model.  I only included a single .amc file from the CMU mocap data base, if you are interested in more mocap data, go to mocap.cs.cmu.edu.  
More generally, if you want to use your own data, pass the data as row-vectors in the Y matrix to gpdmfitFull. 

The code in src/gplvm is curtesy of Neil Lawrence, and src/netlab is simply a copy of the netlab library you can get at http://www.ncrg.aston.ac.uk/netlab/ . 

Cheers,
Jack Wang
