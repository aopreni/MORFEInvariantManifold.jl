function MORFE_integrate_rdyn_frc(analysis_name,zero_amplitude,harmonics_init,param_init,cont_param,
                                  time_integration_length,forward=true,MaxNumPoints=100,minstep=1e-8,maxstep=20.0,ncol=4.0,ntst=40.0,analysis_number=1)
    # check for continuation direction
    if (forward)
      var = 0
    else
      var = 1
    end
    # output folder
    analysis_folder = "./"*analysis_name*"/matcont_automatic"
    analysis_output = "./"*analysis_name*"/frc_"*string(analysis_number)*".mat"
    # 
    nΩ = Int(size(harmonics_init)[1]/2)
    control_parameters = "mu"
    for i = 1:nΩ
      control_parameters *= ",beta"*string(i)
    end
    for i = 1:nΩ
      control_parameters *= ",w"*string(i)
    end
    # initialise normal coordinates and harmonics amplitudes
    X0 = zeros(Float64,size(zero_amplitude)[1]+size(harmonics_init)[1])
    X0[1:size(zero_amplitude)[1]] = zero_amplitude
    X0[size(zero_amplitude)[1]+1:end] = harmonics_init
    # inizialise control parameters
    param = ""
    param *= "mu = "*string(param_init[1])
    for i = 1:nΩ
      param *= ";beta"*string(i) * " = " * string(param_init[i+1])
    end
    for i = 1:nΩ
      param *= ";w"*string(i) * " = " * string(param_init[i+1+nΩ])
    end
    param *= ";"
    println(param)
    #
    mat"""
    warning off
    addpath(genpath('./MatCont7p3')) 
    addpath($analysis_folder)
    ndofs = size($X0);
    ndofs = ndofs(1);
    init 
    X0   = $X0;
    tfin = $time_integration_length;
    """
    eval_string(param)
    println(control_parameters)
    mat"""
    ndofs = size(X0);
    ndofs = ndofs(1);
    hls=feval(@MORFEsystem); 
    options=odeset('RelTol',1e-8);
    """

    param = "[t,y]=ode45(hls{2},[0,tfin],...s
                X0,options,"*control_parameters*");"
    eval_string(param)

    mat"""
    x1 = y(end,:);
    [vl,pk]=findpeaks(y(:,1));
    period=t(pk(end))-t(pk(end-1));
    """

    param = "[t,y] = ode45(hls{2},[0 period],x1,options,"*control_parameters*");"
    eval_string(param)

    mat"""
    active_pars=$cont_param; 
    ncol=$ncol; 
    ntst=$ntst; 
    tolerance=1e-4;
    """
    param = "[x0,v0]=initOrbLC(@MORFEsystem,...
             t,y,...
             ["*control_parameters*"],[active_pars],ntst,ncol,...
             tolerance);"
    eval_string(param)

    mat"""
    opt=contset; 
    opt=contset(opt,'MaxNumPoints',$MaxNumPoints); 
    opt=contset(opt, 'InitStepsize' , 0.1); 
    opt=contset(opt,'MaxStepsize'  , $maxstep); 
    opt=contset(opt,'MinStepsize'  , $minstep); 
    opt=contset(opt,'Backward',0); 
    opt=contset(opt,'FunTolerance', 1e-6);
    opt=contset(opt,'VarTolerance', 1e-6);
    opt=contset(opt,'ActiveParams', active_pars);
    [xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,opt); 
  
    %FRF=[];
    %omega=[];
    %for i=1:size(xlcc,2)
    %   x1=xlcc(1:4:end-2,i);
    %   y1=xlcc(2:4:end-2,i);
    %   omega(i)=(xlcc(end-1,i).^-1)*2*pi;
    %   FRF(i)=max(x1);
    %   FRFv(i)=max(y1);
    %end
    
    %figure(1)
    %hold on
    %plot(omega,FRF)
  
    save($analysis_output,'xlcc','ncol','ntst')
    """
    
    #
    return nothing
    #
  end


function MORFE_integrate_rdyn_backbone(analysis_name,zero_amplitude,time_integration_length,forward=true,MaxNumPoints=100,minstep=1e-8,maxstep=20.0,ncol=4.0,ntst=40.0,analysis_number=1)

  if (forward)
    var = 0
  else
    var = 1
  end

  analysis_folder = "./"*analysis_name*"/matcont_automatic"
  analysis_output = "./"*analysis_name*"/backbone_"*string(analysis_number)*".mat"

  mat"""
  warning off
  addpath(genpath('./MatCont7p3')) 
  addpath($analysis_folder)

  init 

  mu = 0.0;
  ndofs = size($zero_amplitude);
  ndofs = ndofs(1);
  
  X0 = $zero_amplitude;
  
  tfin=$time_integration_length; 
  

  hls=feval(@MORFEsystem); 
  options=odeset('RelTol',1e-8);

  [t,y]=ode45(hls{2},[0,tfin],...s
  X0,options,mu);
  
  x1 = y(end,:);
  
  [vl,pk]=findpeaks(y(:,1));
  period=t(pk(end))-t(pk(end-1));
  [t,y] = ode45(hls{2},[0 period],x1,options,mu);
  
  

  active_pars=[1]; 
  ncol=$ncol; 
  ntst=$ntst; 
  tolerance=1e-4;
  [x0,v0]=initOrbLC(@MORFEsystem,...
  t,y,...
  [mu],active_pars,ntst,ncol,...
  tolerance);
  
  opt=contset; 
  opt=contset(opt,'MaxNumPoints',$MaxNumPoints); 
  opt=contset(opt, 'InitStepsize' , 0.1); 
  opt=contset(opt,'MaxStepsize'  , $maxstep); 
  opt=contset(opt,'MinStepsize'  , $minstep); 
  opt=contset(opt,'Backward',$var); 
  opt=contset(opt,'FunTolerance', 1e-6);
  opt=contset(opt,'VarTolerance', 1e-6);
  [xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,opt); 


  %FRF=[];
  %omega=[];
  %for i=1:size(xlcc,2)
  %    x1=xlcc(1:2:end-2,i);
  %    y1=xlcc(2:2:end-2,i);
  %    omega(i)=(xlcc(end-1,i).^-1)*2*pi;
  %    FRF(i)=max(x1);
  %    FRFv(i)=max(y1);
  %end
  
  %figure(1)
  %hold on
  %plot(omega,FRF)

  save($analysis_output,'xlcc','ncol','ntst')
  
  """

end



function MORFE_compute_backbone_modal(analysis_name,analysis_number=1)
  #
  comb = readdlm(analysis_name*"/comb_a_0/monomials.txt")
  maps = readdlm(analysis_name*"/manifold_a_0/psi.txt")
  #
  po_dir = "./"*analysis_name*"/backbone_"*string(analysis_number)*".mat"
  #
  nc = size(comb)[1]
  ndofs = size(comb)[2]
  neig = size(maps)[2]
  #
  data = matread(po_dir)
  #
  xlcc = data["xlcc"]
  ncol = data["ncol"]
  ntst = data["ntst"]
  #
  npoints = Int(ncol*ntst)
  npo = size(xlcc)[2]
  vars = zeros(Float64,ndofs)
  field = zeros(Float64,neig)
  frf = zeros(npo,neig+1)
  #
  for po = 1:npo
  #
    ω = 2*π/xlcc[end-1,po]
    frf[po,1] = ω
    for p = 1:npoints
      # extract variables
      fill!(field,0.0)
      #
      for d = 1:ndofs
        vars[d] = xlcc[d+(p-1)*ndofs,po]
      end
      #
      for c = 1:nc
        #
        Πi = 1.0
        for d = 1:ndofs 
          Πi *= vars[d]^comb[c,d]
        end
        #
        for d = 1:neig
          field[d] += maps[c,d]*Πi
        end
        #
      end
      #
      for d = 1:neig
        if (abs(field[d])>frf[po,1+d])
          frf[po,1+d] = abs(field[d])
        end
      end
      #
    end
    #
  end
  #
  return frf
  #
end

# reconstruction with only autonomous vector till the paper on the forcing is
# becomes public
function MORFE_compute_frc_modal(analysis_name,Ω_list,analysis_number=1)
  #
  comb = readdlm(analysis_name*"/comb_a_0/monomials.txt")
  maps = readdlm(analysis_name*"/manifold_a_0/psi.txt")
  #
  po_dir = "./"*analysis_name*"/frc_"*string(analysis_number)*".mat"
  #
  nc = size(comb)[1]
  ndofs = size(comb)[2]
  nΩ = size(Ω_list)[1]
  neig = size(maps)[2]
  #
  data = matread(po_dir)
  #
  xlcc = data["xlcc"]
  ncol = data["ncol"]
  ntst = data["ntst"]
  #
  npoints = Int(ncol*ntst)
  npo = size(xlcc)[2]
  vars = zeros(Float64,ndofs)
  field = zeros(Float64,neig)
  frf = zeros(npo,neig+1)
  #
  for po = 1:npo
  #
    ω = xlcc[end,po]
    frf[po,1] = ω
    for p = 1:npoints
      # extract variables
      fill!(field,0.0)
      #
      for d = 1:ndofs
        vars[d] = xlcc[d+(p-1)*(ndofs+nΩ),po]
      end
      #
      for c = 1:nc
        #
        Πi = 1.0
        for d = 1:ndofs 
          Πi *= vars[d]^comb[c,d]
        end
        #
        for d = 1:neig
          field[d] += maps[c,d]*Πi
        end
        #
      end
      #
      for d = 1:neig
        if (abs(field[d])>frf[po,1+d])
          frf[po,1+d] = abs(field[d])
        end
      end
      #
    end
    #
  end
  #
  return frf
  #
end
