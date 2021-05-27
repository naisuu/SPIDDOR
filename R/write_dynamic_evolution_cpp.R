write.dynamic_evolution_cpp=function(arguments,fun_header,nodes.names)
{
  #write("","dyn_evolution.cpp")
 
  Head<-paste("
// [[Rcpp::export]]
NumericVector time_evolution_f(const int& ts",
              arguments,",std::vector<std::string> Knockouts,std::vector<std::string> Over_expr,
              std::vector<std::string> Over_expr_AA, std::vector<std::string> Mutations,
              List MCAB_times,List Tx_times, List MUT_times, List Polym,
              IntegerVector Initial_cond, const bool asynchronous){

  std::string node_i;
  std::vector<int> samples;")
  
  nodes<-paste("\tstd::string nodes_names[] = {",paste('"',paste(nodes.names,collapse='","'),'"',sep=""),"};",sep="") #Nombre de los nodos entre comillas
  
  dyn1<-"
  const int n_nodes = sizeof(nodes_names) / sizeof(nodes_names[0]);
  // Polymorphism constant when applying polymorphism for specific time steps only
  const float polym_val = 0.1; // TODO: Should I make this user input? Eventually can also be applied with different values for different nodes, e.g. Mutations=c("A", "B"), Vals=c(0.1, 0.5)
  //Pattern and update creation:
  int *pattern = new int[n_nodes*(ts + 1)]();
  int *update =new int[n_nodes*(ts + 1)]();
          
  //Initial conditions:
  for (int it = 0; it < Initial_cond.size(); it++){ pattern[Initial_cond[it]*(ts + 1)]=(1);}
  for (int i = 0; i < n_nodes; i++) { update[i*(ts)+i] = ( 1 ); }\n"
  Polymorphism<-c()
  for(i in nodes.names)
    Polymorphism<-c(Polymorphism,paste("\tdouble P_",i,'=as<double>(Polym["',i,'"]);',sep="") ) 
  
  Polyms_to_update <-c()
  for(i in nodes.names)
    Polyms_to_update <-c(Polyms_to_update, paste("\tP_",i,'=as<double>(Polym["',i,'"]);',sep=""))
  
  dyn2<-"\n  
  //Iterate:
  for (int j = 1; j <= ts; j++) {
    samples = myrandom(n_nodes);
    for (int i = 0; i <n_nodes; i++) {
      node_i = nodes_names[samples.at(i)];
              
      if (std::find(Over_expr_AA.begin(), Over_expr_AA.end(), node_i) != Over_expr_AA.end()) {
        int node = std::distance(nodes_names, std::find(nodes_names, nodes_names + (n_nodes - 1), node_i));
        if (pattern[node*(ts + 1) + (j - 1)] == 1) {
          pattern[node*(ts + 1) + j] = 1;
          continue;
        }
      }
      
      // Check if mutations (polymorphisms) are passed as arguments
      if (std::find(Mutations.begin(), Mutations.end(), node_i) != Mutations.end()) {
        // If mutations were requested but the timesteps were not, this mutation applies indefinitely
        // This is kind of redundant as this feature already existed, however this can serve as syntax sugar
        if (MUT_times.size() == 0) {
          Polym[node_i] = polym_val; // Set the polymorphism of that node to our preset polymorphism constant
        } else {
          // if a time step range is given, only set the polymorphism constant for those time steps.
          int pos = std::find(Mutations.begin(), Mutations.end(), node_i) - Mutations.begin();
          std::vector<int>MUT_time = MUT_times[pos];
          
          if (std::find(MUT_time.begin(), MUT_time.end(), j) != MUT_time.end()) {
            Polym[node_i] = polym_val;
          } else {
            // Once the specified time steps have been completed, set the polymorphism back to default (1)
            Polym[node_i] = 1; // TODO: Allow users to combine this with preset polymorphism, i.e. store original initial conditions separately.
          }
        }
      }
      
      
      if (std::find(Knockouts.begin(), Knockouts.end(), node_i) != Knockouts.end()){// && MCAB_times.size()==0){
        int node = std::distance(nodes_names, std::find(nodes_names, nodes_names + (n_nodes - 1), node_i));
        if(MCAB_times.size()==0){
          pattern[node*(ts + 1) + j] = 0;
          update[samples.at(i)*(ts+1) + j] = 1;
          continue;
        }else{
          int pos = std::find(Knockouts.begin(), Knockouts.end(), node_i) - Knockouts.begin();
          std::vector<int>KO_time=MCAB_times[pos];
          if(std::find(KO_time.begin(), KO_time.end(), j) != KO_time.end()){
            pattern[node*(ts + 1) + j] = 0;
            update[samples.at(i)*(ts+1) + j] = 1;
            continue;
          }
        }
      }
          
      if (std::find(Over_expr.begin(), Over_expr.end(), node_i) != Over_expr.end()){
        int node = std::distance(nodes_names, std::find(nodes_names, nodes_names + (n_nodes - 1), node_i));
        if(Tx_times.size()==0){
          pattern[node*(ts + 1) + j] = 1;
          update[samples.at(i)*(ts+1) + j] = 1;
          continue;
        }else{
          int pos = std::find(Over_expr.begin(), Over_expr.end(), node_i) - Over_expr.begin();
          std::vector<int>OE_time=Tx_times[pos];
          if(std::find(OE_time.begin(), OE_time.end(), j) != OE_time.end()){
            pattern[node*(ts + 1) + j] = 1;
            update[samples.at(i)*(ts+1) + j] = 1;
            continue;
          }
        }
      }
  "      
  
  Pie_dyn<-"
      update[samples.at(i)*(ts+1) + j] = 1;
    }
  } 
  Rcpp::NumericVector P(pattern, pattern +(n_nodes*(ts + 1)));
  delete[] update;
  delete[] pattern;
  return(P);
}"
  All<-c(Head,nodes,dyn1,Polymorphism,dyn2,Polyms_to_update,fun_header,Pie_dyn)
  write(All,"Boolean_func_C.cpp",append=TRUE)
}

