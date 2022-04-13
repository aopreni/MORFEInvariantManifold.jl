
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
