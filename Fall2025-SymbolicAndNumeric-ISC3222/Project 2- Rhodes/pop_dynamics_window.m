function pop_dynamics_tabbed_colorblind
% Full runnable population dynamics simulator with tabs, food regrowth control,
% role-specific Assessor behavior, and colorblind-friendly colors

%% ---------------- FIGURE ----------------
f = figure('Name','Population Dynamics Simulator','Position',[100 100 1200 600]);

% Create tab group
tg = uitabgroup(f,'Position',[0 0 1 1]);

% Description tab
tabDesc = uitab(tg,'Title','Simulation Description');
descStr = ['This simulation models pigeon population dynamics using four roles:' newline ...
    '1) Hawks (HH): aggressive, always fight.' newline ...
    '2) Doves (DD): passive, avoid fights.' newline ...
    '3) Retaliators (HD): fight only if attacked.' newline ...
    '4) Assessors (AS): fight if stronger than opponent.' newline newline ...
    'Genetics:' newline ...
    '- Each bird has 2 loci: aggression and perception.' newline ...
    '- Offspring inherit one allele per locus from each parent (Mendelian inheritance).' newline ...
    '- Genetics influence role behavior if "Enable Genetics" is checked.' newline newline ...
    'Reproduction and death are based on fitness.' newline ...
    'Food can be added to sandbox by clicking.'];
uicontrol('Parent',tabDesc,'Style','edit','Max',2,'Enable','inactive', ...
    'Units','normalized','Position',[0.05 0.05 0.9 0.9],'String',descStr,'FontSize',11);

% Simulation tab
tabSim = uitab(tg,'Title','Run Simulation');

%% ---------------- CONTROLS PANEL ----------------
panel = uipanel('Parent',tabSim,'Title','Population Controls','FontSize',11,'Position',[0.02 0.05 0.22 0.9]);

y = 0.93; dy = 0.06;

uicontrol(panel,'Style','text','String','Initial Populations','FontWeight','bold','FontSize',10,'Units','normalized','Position',[0.05 y 0.9 0.05]);

y = y - dy;
uicontrol(panel,'Style','text','String','Hawks (HH):','Units','normalized','Position',[0.05 y 0.45 0.05]);
hawkBox = uicontrol(panel,'Style','edit','String','10','Units','normalized','Position',[0.55 y 0.35 0.05]);

y = y - dy;
uicontrol(panel,'Style','text','String','Doves (DD):','Units','normalized','Position',[0.05 y 0.45 0.05]);
doveBox = uicontrol(panel,'Style','edit','String','10','Units','normalized','Position',[0.55 y 0.35 0.05]);

y = y - dy;
uicontrol(panel,'Style','text','String','Retaliators (HD):','Units','normalized','Position',[0.05 y 0.45 0.05]);
retBox = uicontrol(panel,'Style','edit','String','10','Units','normalized','Position',[0.55 y 0.35 0.05]);

y = y - dy;
uicontrol(panel,'Style','text','String','Assessors (AS):','Units','normalized','Position',[0.05 y 0.45 0.05]);
assBox = uicontrol(panel,'Style','edit','String','10','Units','normalized','Position',[0.55 y 0.35 0.05]);

y = y - dy*1.3;
genCheck = uicontrol(panel,'Style','checkbox','String','Enable Genetics','Value',1,'Units','normalized','Position',[0.05 y 0.9 0.05]);

y = y - dy*1.1;
uicontrol(panel,'Style','text','String','Prob. of Reproduction (0–1):','Units','normalized','Position',[0.05 y 0.9 0.05]);
y = y - dy;
reproBox = uicontrol(panel,'Style','edit','String','0.2','Units','normalized','Position',[0.05 y 0.4 0.05]);

y = y - dy*1.2;
uicontrol(panel,'Style','text','String','Offspring Starting Fitness:','Units','normalized','Position',[0.05 y 0.9 0.05]);
y = y - dy;
offspringBox = uicontrol(panel,'Style','edit','String','20','Units','normalized','Position',[0.05 y 0.4 0.05]);

y = y - dy*1.2;
uicontrol(panel,'Style','text','String','Food regrowth per timestep:','Units','normalized','Position',[0.05 y 0.9 0.05]);
y = y - dy;
foodRegrowBox = uicontrol(panel,'Style','edit','String','5','Units','normalized','Position',[0.05 y 0.4 0.05]);

y = y - dy*1.5;
uicontrol(panel,'Style','pushbutton','String','Start Simulation','Units','normalized','FontSize',10,...
    'Position',[0.1 y 0.8 0.07],'Callback',@runSimulation);

%% ---------------- PLOTS ----------------
axPop = axes('Parent',tabSim,'Position',[0.30 0.1 0.30 0.8]);
title(axPop,'Population Over Time'); xlabel(axPop,'Time'); ylabel(axPop,'Count'); hold(axPop,'on');

axSand = axes('Parent',tabSim,'Position',[0.64 0.1 0.32 0.8]);
title(axSand,'Sandbox (Birds + Food)'); axis(axSand,[0 1 0 1]); axis(axSand,'equal'); hold(axSand,'on');
axSand.XTick=[]; axSand.YTick=[];

food=[]; 
set(f,'WindowButtonDownFcn',@clickAddFood);

%% ---------------- CALLBACKS ----------------
function clickAddFood(~,~)
    clickFig = get(f,'CurrentPoint');
    figPos = get(f,'Position');
    clickNorm = [clickFig(1)/figPos(3), clickFig(2)/figPos(4)];
    axPos = get(axSand,'Position');
    insideX = clickNorm(1)>=axPos(1) && clickNorm(1)<=axPos(1)+axPos(3);
    insideY = clickNorm(2)>=axPos(2) && clickNorm(2)<=axPos(2)+axPos(4);
    if ~(insideX && insideY), return; end
    axX=(clickNorm(1)-axPos(1))/axPos(3);
    axY=(clickNorm(2)-axPos(2))/axPos(4);
    nFood=25;jitter=0.03;
    newFood=[axX+jitter*(rand(nFood,1)-0.5), axY+jitter*(rand(nFood,1)-0.5)];
    newFood(newFood<0)=0; newFood(newFood>1)=1;
    food=[food; newFood];
end

%% ---------------- SIMULATION ----------------
function runSimulation(~,~)
    % Read inputs
    pop_HH = str2double(hawkBox.String);
    pop_DD = str2double(doveBox.String);
    pop_HD = str2double(retBox.String);
    pop_AS = str2double(assBox.String);
    P_repro = str2double(reproBox.String);
    startFitness = str2double(offspringBox.String);
    foodRegrow = str2double(foodRegrowBox.String);
    useGenetics = genCheck.Value;

    % Initialize individuals
    genotypes = [repmat("HH",1,pop_HH), repmat("DD",1,pop_DD), repmat("HD",1,pop_HD), repmat("AS",1,pop_AS)];
    fitness = startFitness*ones(1,numel(genotypes));
    positions = rand(numel(genotypes),2);
    velocities = zeros(numel(genotypes),2);

    % Initialize food
    targetFood = max(50,3*numel(genotypes));
    food = rand(targetFood,2);

    baseSpeed = 0.02;
    collisionRadius = 0.05;
    timesteps = 300;
    history = zeros(timesteps,4);

    % Genetics: 2 loci per bird (aggression, perception)
    loci = randi([0 1],numel(genotypes),2); % 0/1 alleles

    % Role colors (colorblind-friendly)
    roleColors = containers.Map({'HH','DD','HD','AS'},{[1 0.5 0],[0 0 1],[0.5 0 0.5],[0 1 1]});

    % Main loop
    for t = 1:timesteps
        if isempty(genotypes), break; end
        N = numel(genotypes);

        % Movement: chase nearest food
        for i=1:N
            if isempty(food)
                dir = rand(1,2)-0.5;
            else
                dif = food - positions(i,:);
                d2 = sum(dif.^2,2);
                [~,idx] = min(d2);
                target = dif(idx,:);
                dir = target/norm(target);
            end
            velocities(i,:) = dir*baseSpeed;
        end

        positions = positions + velocities;
        % Bounce edges
        for d=1:2
            over = positions(:,d)>1;
            under = positions(:,d)<0;
            positions(over,d)=1; velocities(over,d)=-abs(velocities(over,d));
            positions(under,d)=0; velocities(under,d)=abs(velocities(under,d));
        end

        % Collision detection
        collPair = [];
        for i=1:N-1
            for j=i+1:N
                if norm(positions(i,:)-positions(j,:))<collisionRadius
                    collPair = [i j]; break;
                end
            end
            if ~isempty(collPair), break; end
        end

        if ~isempty(collPair)
            i1=collPair(1); i2=collPair(2);
            s1=genotypes(i1); s2=genotypes(i2);

            % Simple payoff matrix
            V=50; C=100;
            M = [ (V-C)/2 V (V-C)/2 V;
                  0      V/2 V/2  V;
                  (V-C)/2 V  V/2  V;
                  (V-C)/2 V  V/2  V];

            % --- Role-specific Assessor behavior ---
            if s1=="AS" && fitness(i1)<fitness(i2)
                payoff1 = 0; payoff2 = 0;
            else
                payoff1 = M(strategy_index(s1),strategy_index(s2));
                payoff2 = M(strategy_index(s2),strategy_index(s1));
            end
            if s2=="AS" && fitness(i2)<fitness(i1)
                payoff1 = 0; payoff2 = 0;
            end

            fitness(i1)=fitness(i1)+payoff1;
            fitness(i2)=fitness(i2)+payoff2;

            % Reproduction
            if rand<P_repro
                if useGenetics
                    baby = mendelian(loci(i1,:),loci(i2,:),s1,s2);
                else
                    baby.geno=s1;
                    baby.loci=loci(i1,:);
                end
                genotypes(end+1)=baby.geno;
                loci(end+1,:)=baby.loci;
                fitness(end+1)=startFitness;
                positions(end+1,:)=(positions(i1,:)+positions(i2,:))/2 + 0.01*(rand(1,2)-0.5);
                velocities(end+1,:)=[0 0];
            end

            % Remove nearest food
            if ~isempty(food)
                midpt=(positions(i1,:)+positions(i2,:))/2;
                dif=food-midpt;
                d2=sum(dif.^2,2);
                [~,idx]=min(d2);
                food(idx,:)=[];
            end
        end

        % Death
        alive = fitness>0;
        genotypes = genotypes(alive);
        fitness = fitness(alive);
        positions = positions(alive,:);
        velocities = velocities(alive,:);
        loci = loci(alive,:);

        % Food regrowth
        currFood=size(food,1);
        if currFood<targetFood
            toAdd=min(foodRegrow,targetFood-currFood);
            food=[food; rand(toAdd,2)];
        end

        % Update population plot
        countH=sum(genotypes=="HH"); countD=sum(genotypes=="DD");
        countR=sum(genotypes=="HD"); countA=sum(genotypes=="AS");
        history(t,:)=[countH countD countR countA];
        cla(axPop)
        plot(axPop,history(1:t,1),'Color',roleColors('HH'),'LineWidth',2); hold(axPop,'on')
        plot(axPop,history(1:t,2),'Color',roleColors('DD'),'LineWidth',2)
        plot(axPop,history(1:t,3),'Color',roleColors('HD'),'LineWidth',2)
        plot(axPop,history(1:t,4),'Color',roleColors('AS'),'LineWidth',2)
        legend(axPop,{'HH','DD','HD','AS'},'Location','northwest');

        % Update sandbox
        cla(axSand); hold(axSand,'on'); axis(axSand,[0 1 0 1]); axis(axSand,'equal')
        if ~isempty(food), scatter(axSand,food(:,1),food(:,2),10,'k','filled'); end
        if ~isempty(genotypes)
            cols=zeros(numel(genotypes),3);
            for k=1:numel(genotypes)
                cols(k,:)=roleColors(genotypes(k));
            end
            scatter(axSand,positions(:,1),positions(:,2),80,cols,'filled');
        end
        drawnow
    end
end

%% ---------------- HELPER FUNCTIONS ----------------
function idx = strategy_index(s)
    switch s
        case "HH", idx=1;
        case "DD", idx=2;
        case "HD", idx=3;
        case "AS", idx=4;
    end
end

function baby = mendelian(loci1,loci2,g1,g2)
    newLoci = [randsample([loci1(1) loci2(1)],1) randsample([loci1(2) loci2(2)],1)];
    if rand<0.5, newG=g1; else newG=g2; end
    baby.geno=newG;
    baby.loci=newLoci;
end

end




