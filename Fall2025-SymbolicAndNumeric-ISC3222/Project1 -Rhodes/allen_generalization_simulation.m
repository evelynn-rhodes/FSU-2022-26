function allen_generalization_viewer()
% ========================================================================
% Allen Generalization Viewer
% Mode 1 -> Behavioral + neural simulation using generalization G
% Mode 2 -> Test a real/specified neural response for generalization
% ========================================================================

clc; disp("Allen Neural Generalization System");
disp("1 = Simulation mode (neural manifold + G)");
disp("2 = Single-trial neural generalization test");

mode = input("Select mode (1/2): ");
if isempty(mode), mode = 1; end

switch mode
    case 1
        run_simulation();
    case 2
        run_neuron_test();
    otherwise
        error("Invalid mode.");
end
end


%% ========================================================================
%% ============= MODE 1 — BEHAVIORAL + NEURAL SIMULATION ==================
%% ========================================================================

function run_simulation()
% Uses Resp_fixed to simulate neural activity & project into PCA space

%% ---- Check data ----
req = {'Resp_fixed','cell_ids','ori_expected','sf_expected','exp_id'};
for r=req
    assert(evalin("base",sprintf("exist('%s','var')",r{1}))==1, ...
        "Missing %s in workspace",r{1});
end

Resp_fixed   = evalin("base","Resp_fixed");   % [6×5×cells]
cell_ids     = evalin("base","cell_ids");
ori_expected = evalin("base","ori_expected");
sf_expected  = evalin("base","sf_expected");
exp_id       = evalin("base","exp_id");
[nOri,nSF,nCell] = size(Resp_fixed);

%% ---- ORI mapping ----
target_oris = [0 30 60 90 120 150];
ori_map = zeros(1,6);
for i=1:6
    oi = find(ori_expected==target_oris(i),1);
    assert(~isempty(oi),"Orientation %d missing",target_oris(i));
    ori_map(i)=oi;
end

%% ---- Build templates ----
templates_raw = zeros(nCell,6);
for i=1:6
    T = squeeze(Resp_fixed(ori_map(i),:,:)).';    % cell × SF
    templates_raw(:,i) = mean(T,2,'omitnan');
end

mu = mean(templates_raw,2);
sd = std(templates_raw,0,2); sd(sd==0)=1;
templates_z = (templates_raw - mu) ./ sd;

Lidx = [1 2 3]; Ridx = [4 5 6];
protoL = mean(templates_z(:,Lidx),2);
protoR = mean(templates_z(:,Ridx),2);

%% ---- Input generalization coefficient ----
G = input("Enter generalization coefficient (0–1, default 0.7): ");
if isempty(G), G = 0.7; end
G = max(0,min(1,G));
noise_sigma = 1.0;
amp_jitter  = 0.15;
softmax_temp = 1.8;

%% ---- Simulation sequence ----
nTrials = 25;
ori_seq = target_oris(randi(6,nTrials,1));
truth = containers.Map(num2cell(target_oris),[repmat({'L'},1,3),repmat({'R'},1,3)]);
correct_count = 0;

%% ---- Build PCA space (30 stimuli × cells) ----
stimMat = reshape(permute(Resp_fixed,[1 2 3]),[nOri*nSF,nCell]);
col_mu = mean(stimMat,1);
col_sd = std(stimMat,0,1); col_sd(col_sd==0)=1;
X = (stimMat-col_mu)./col_sd;
[coeff,score,~,~,~,mu_pca] = pca(X);
stimPC = score(:,1:2);

project = @(x) (((x(:)'-col_mu)./col_sd - mu_pca)*coeff(:,1:2));

sf_colors = lines(nSF);
ori_idx = repelem((1:nOri)',nSF);
sf_idx  = repmat((1:nSF)',nOri,1);

%% ---- UI ----
fig = figure("Color","w","Position",[100 80 1200 720]);
tiledlayout(fig,2,3,"TileSpacing","compact","Padding","compact");

btn_prev = uicontrol("style","pushbutton","string","⟵ Prev", ...
    "units","normalized","position",[0.01 0.95 0.08 0.04], ...
    "callback",@(~,~) step(-1));

btn_next = uicontrol("style","pushbutton","string","Next ⟶", ...
    "units","normalized","position",[0.10 0.95 0.08 0.04], ...
    "callback",@(~,~) step(+1));

btn_play = uicontrol("style","togglebutton","string","▶ Play", ...
    "units","normalized","position",[0.19 0.95 0.08 0.04], ...
    "callback",@(~,~) playLoop);

trial = 1; playing = false;
updateFrame();

%% ---- Frame functions ----
function playLoop(~,~)
    playing = get(btn_play,'value');
    if playing, set(btn_play,'string',"⏸ Pause"); end
    while playing
        pause(0.7)
        trial = min(trial+1,nTrials);
        updateFrame();
        playing = get(btn_play,'value');
    end
    set(btn_play,'value',0,'string',"▶ Play");
end

function step(d)
    trial = max(1,min(nTrials,trial+d));
    updateFrame();
end

function updateFrame()
    tiledlayout(fig,2,3);

    ori = ori_seq(trial);
    trueSide = truth(ori);
    oi = find(target_oris==ori);

    base = templates_raw(:,oi);
    noise = noise_sigma*randn(nCell,1);
    jitter = 1+amp_jitter*randn;
    x_raw = (G*jitter)*base + (1-G)*noise;

    x = (G*jitter)*templates_z(:,oi) + (1-G)*randn(nCell,1);

    simL = dot(x,protoL)/(norm(x)*norm(protoL));
    simR = dot(x,protoR)/(norm(x)*norm(protoR));
    simL = G*simL; simR = G*simR;
    pR = exp(simR/softmax_temp)/(exp(simL/softmax_temp)+exp(simR/softmax_temp));
    pred = "R"; if rand>=pR, pred="L"; end
    correct = (pred == trueSide);

    if correct, correct_count = correct_count + 1; end
    pt = project(x_raw);

    %% Stim
    nexttile(1); cla; axis off;
    title(sprintf("Stim %d°",ori),"FontWeight","bold");
    drawOri(ori);

    %% Text
    nexttile(2); cla; axis off;
    text(0,0.9,sprintf("Trial %d/%d",trial,nTrials),"FontSize",14,"FontWeight","bold");
    text(0,0.75,sprintf("G = %.2f",G),"FontSize",12);
    text(0,0.62,sprintf("Truth = %s | Pred = %s",trueSide,pred),"FontSize",12);
    text(0,0.50,sprintf("p(R)=%.2f",pR),"FontSize",12);
    text(0,0.35,sprintf("Correct: %d/%d",correct_count,trial),"FontWeight","bold");

    %% PCA
    nexttile(3); cla; hold on;
    alpha=.35; lw=2; ms=7;
    for sf = 1:nSF
        idx = (sf_idx==sf);
        pts = stimPC(idx,:);
        plot(pts(:,1),pts(:,2),"--","Color",[sf_colors(sf,:) alpha],"LineWidth",lw);
        scatter(pts(:,1),pts(:,2),ms,"filled", ...
            "MarkerFaceColor",sf_colors(sf,:), ...
            "MarkerEdgeColor","k","MarkerFaceAlpha",0.6);
    end
    scatter(pt(1),pt(2),120,'p',"filled","MarkerFaceColor","y","MarkerEdgeColor","k");
    xlabel("PC1"); ylabel("PC2"); title("Neural state in PCA"); axis equal; grid on;

    %% Rule
    nexttile(6); cla; axis off;
    text(0,0.9,"Training Rule","FontWeight","bold");
    text(0,0.7,"Left: 0,30,60°");
    text(0,0.55,"Right: 90,120,150°");
    text(0,0.3,sprintf("Cells: %d",nCell));
    text(0,0.15,sprintf("Exp %d",exp_id));
end

end % simulation




%% ========================================================================
%% ============== MODE 2 — REAL NEURAL INPUT CHECK ========================
%% ========================================================================

function run_neuron_test()
% Tests whether a user-provided or synthetic ΔF/F trial
% matches (generalizes to) the neural manifold.

%% ---- check data ----
req={'Resp_fixed','ori_expected','sf_expected','cell_ids','exp_id'};
for r=req
    assert(evalin("base",sprintf("exist('%s','var')",r{1}))==1, ...
        "Missing %s",r{1});
end

Resp_fixed   = evalin("base","Resp_fixed");
ori_expected = evalin("base","ori_expected");
sf_expected  = evalin("base","sf_expected");
cell_ids     = evalin("base","cell_ids");
exp_id       = evalin("base","exp_id");
[nOri,nSF,nCell] = size(Resp_fixed);

%% ---- Get user inputs ----
ori_in = input("Orientation (0–150°): ");
if isempty(ori_in), ori_in=ori_expected(randi(nOri)); end

sf_in = input("Spatial frequency (0.02–0.32): ");
if isempty(sf_in), sf_in=sf_expected(randi(nSF)); end

cell_in = input("Cell ID (blank = use all cells): ");
if ~isempty(cell_in)
    cidx = find(cell_ids==cell_in,1);
    assert(~isempty(cidx),"Cell ID not found");
else
    cidx=[];
end

trace = input("Enter ΔF/F trace (blank = auto): ");
if isempty(trace)
    [~,oi]=min(abs(ori_expected-ori_in));
    [~,si]=min(abs(sf_expected-sf_in));
    base = squeeze(Resp_fixed(oi,si,:));
    trace = mean(base)*0.01 + randn(1,600)*0.01;
end
trace = trace(:)';

%% ---- compute ΔF/F response like Allen ----
b = mean(trace(1:round(0.2*numel(trace))));
resp_trial = mean(trace) - b;

if ~isempty(cidx)
    x = resp_trial;
else
    Resp_avg=squeeze(mean(Resp_fixed,2));
    [~,oi]=min(abs(ori_expected-ori_in));
    base=Resp_avg(oi,:)';
    x = base + resp_trial + 0.05*randn(nCell,1);
end

%% ---- load PCA space ----
stimMat = reshape(permute(Resp_fixed,[1 2 3]),[nOri*nSF,nCell]);
mu = mean(stimMat,1); sd = std(stimMat,0,1); sd(sd==0)=1;
X=(stimMat-mu)./sd;
[coeff,score,~,~,~,mu_pca]=pca(X);
stimPC = score(:,1:2);
project=@(v)(((v(:)'-mu)./sd - mu_pca)*coeff(:,1:2));
pt = project(x);

%% ---- distance & decision ----
D = vecnorm(stimPC-pt,2,2);
dmin = min(D);
th = input("Generalization threshold (default 0.65): ");
if isempty(th), th=0.65; end
generalized = dmin < th;

%% ---- output ----
fprintf("\nOrientation = %.1f°, SF = %.3f\n",ori_in,sf_in);
fprintf("Nearest manifold distance = %.4f\n",dmin);
fprintf("Threshold = %.4f\n",th);
if generalized
    fprintf("✅ Neural response GENERALIZED\n");
else
    fprintf("❌ Neural response DID NOT GENERALIZE\n");
end

%% ---- plot ----
figure("Color","w","Position",[100 100 600 450]); hold on
sf_colors = lines(nSF);
ori_idx = repelem((1:nOri)',nSF);
sf_idx  = repmat((1:nSF)',nOri,1);

alpha=0.35; lw=2; ms=7;
for sf=1:nSF
    idx=(sf_idx==sf);
    pts=stimPC(idx,:);
    plot(pts(:,1),pts(:,2),"--","Color",[sf_colors(sf,:) alpha],"LineWidth",lw);
    scatter(pts(:,1),pts(:,2),ms,"filled", ...
        "MarkerFaceColor",sf_colors(sf,:), ...
        "MarkerEdgeColor","k","MarkerFaceAlpha",0.6);
end

scatter(pt(1),pt(2),160,'p',"filled","MarkerFaceColor","y","MarkerEdgeColor","k","LineWidth",1.5);
title("Neural trial projection in PCA manifold");
xlabel("PC1"); ylabel("PC2"); grid on; axis equal;

end % neuron test



%% ========================================================================
%% helper: draw oriented grating
function drawOri(ori)
[X,Y]=meshgrid(linspace(0,1,200));
th=deg2rad(ori);
img=sin(2*pi*6*(X*cos(th)+Y*sin(th)));
imagesc(img); colormap gray; axis equal off tight;
end

