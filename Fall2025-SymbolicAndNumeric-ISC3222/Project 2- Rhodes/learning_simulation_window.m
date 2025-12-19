function learning_simulations_window
%% LEARNING SIMULATIONS FOR PIGEONS
% Single Learner, Multi-Agent, and Best Role Evaluation
% All 4 roles: Hawk, Dove, Retaliator, Assessor

%% ---------------- MAIN FIGURE -----------------
fig = uifigure('Name','Learning Simulations','Position',[100 100 1000 600]);

tg = uitabgroup(fig,'Position',[0 0 fig.Position(3) fig.Position(4)]);

% ------------------- LEFT TAB -------------------
tabText = uitab(tg,'Title','Info');
% Use uitextarea if you want scrollable text
uitextarea(tabText,'Value', [
        "Welcome to the Learning Simulation!"
        "-----------------------------------"
        "This tab contains plain text explaining the simulation, roles, and learning concepts."
        "Learning allows for the pigeons to switch between roles to find the most suitable role for them."
        "In the learning simulations, everything really comes down to two main variables. The first is alpha (α), which is the learning rate. Alpha controls how quickly the pigeon updates what it has learned — a high alpha means it changes its mind fast after each interaction, while a low alpha means it updates slowly and stays more steady. The second key variable is tau (τ), which controls how random or exploratory the pigeon’s choices are. A high tau makes the pigeon try different strategies more often, and a low tau makes it consistently pick whatever it currently thinks is best. All of these learning models revolve around something called Q-values, which are basically the pigeon’s “beliefs” or expectations about how good each strategy is. Each time the pigeon gets a reward, it updates that belief using the core equation Q(new) = Q(old) + α × (reward − Q(old)), meaning the pigeon adjusts its expectations based on how surprising the outcome was. Then, to decide which action to take, it uses a softmax equation that turns those Q-values into probabilities based on tau. These two variables, combined with these simple equations, are the cornerstone of how the learning simulations work. "
        "In the individual simulation, a single pigeon uses these Q-values to learn which of the four strategies works best by repeatedly playing against randomly chosen opponents. Each round, it picks a strategy using the softmax rule, gets a reward from the matrix M, and then updates the Q-value for the strategy it just used with the learning rule described above. Over time, strategies that earn better rewards naturally get higher Q-values, so the pigeon chooses them more often. The plots simply show how its expectations and choice probabilities change as it learns."
        "In the multi-agent simulation, instead of one pigeon learning on its own, you now have a whole group of pigeons (N agents) learning at the same time. Each agent keeps its own set of Q-values and, on every round, each one chooses a strategy using the same softmax rule based on tau. After selecting a strategy, the agent is matched with a randomly chosen opponent strategy and receives a reward from the matrix M, then updates its Q-value for the chosen action using the same learning rule as before. Because all agents are learning in parallel, the simulation tracks how the average Q-values across the population change over time and how often each strategy is being chosen in the group. The result is a picture of how learning spreads through a population, showing which strategies rise or fall in popularity as the agents collectively gain experience."
        "In this simulation, each of the N learners tracks a Q-value for each of the four strategies, and these Q-values represent how good each strategy has seemed based on past experience. On every trial, each learner uses a softmax choice rule, where their current Q-values are divided by the temperature tau, exponentiated, and normalized so that higher Q’s are slightly (or strongly, if tau is small) more likely to be chosen; this produces a probability distribution the agent samples an action from. After choosing a strategy a and facing a randomly selected opponent strategy, the agent receives a reward r from the payoff matrix M, and then updates only the chosen Q-value using the learning rule Q(a) ← Q(a) + α(r − Q(a)), where alpha controls how quickly the agent learns from new outcomes. Over many trials, this process causes the Q-values—and therefore the strategy choices—to shift toward whichever strategies consistently produce better rewards."
        "The first two simulations also plots the probability of choosing a role (policy probabilities), just because a pigeon belives that one role will do the best, depending on the value of tau it may decide to choose a different role instead to explore if that role is better than the role with the highest Q-value."
    ], 'Editable',false, 'Position',[20 20 760 560]);

%% ---------------- TAB 1: SINGLE LEARNER -----------------
tab1 = uitab(tg,'Title','Single Learner');
uilabel(tab1,'Text',['A single pigeon learns via Q-learning which role to pick. ',...
    'Plots show Q-values and policy probabilities over time.'],...
    'Position',[20 500 950 50],'FontSize',14);

ctrl1 = uipanel(tab1,'Title','Controls','Position',[10 10 250 480]);
% Controls
uilabel(ctrl1,'Text','Learning Rate α:','Position',[10 430 120 20]);
lrBox1 = uieditfield(ctrl1,'numeric','Value',0.2,'Limits',[0 1],'Position',[140 430 80 20]);
uilabel(ctrl1,'Text','Softmax τ:','Position',[10 390 120 20]);
tauBox1 = uieditfield(ctrl1,'numeric','Value',1,'Position',[140 390 80 20]);
uilabel(ctrl1,'Text','Number of Trials:','Position',[10 350 120 20]);
trialBox1 = uieditfield(ctrl1,'numeric','Value',300,'Position',[140 350 80 20]);

uibutton(ctrl1,'Text','Run Simulation','Position',[50 300 150 40],...
    'ButtonPushedFcn',@runSingleLearner);

% Axes
axQ1 = uiaxes(tab1,'Position',[300 320 680 200]);
title(axQ1,'Q-values'); xlabel(axQ1,'Trial'); ylabel(axQ1,'Q'); hold(axQ1,'on'); grid(axQ1,'on');
axP1 = uiaxes(tab1,'Position',[300 50 680 200]);
title(axP1,'Policy Probabilities'); xlabel(axP1,'Trial'); ylabel(axP1,'P(a)'); hold(axP1,'on'); grid(axP1,'on');

%% ---------------- TAB 2: MULTI-AGENT LEARNING -----------------
tab2 = uitab(tg,'Title','Multi-Agent Learning');
uilabel(tab2,'Text',['Multiple pigeons (agents) learn simultaneously. ',...
    'Track how strategies evolve in a population over trials.'],...
    'Position',[20 500 950 50],'FontSize',14);

ctrl2 = uipanel(tab2,'Title','Controls','Position',[10 10 250 480]);
uilabel(ctrl2,'Text','Learning Rate α:','Position',[10 430 120 20]);
lrBox2 = uieditfield(ctrl2,'numeric','Value',0.2,'Limits',[0 1],'Position',[140 430 80 20]);
uilabel(ctrl2,'Text','Softmax τ:','Position',[10 390 120 20]);
tauBox2 = uieditfield(ctrl2,'numeric','Value',1,'Position',[140 390 80 20]);
uilabel(ctrl2,'Text','Number of Agents:','Position',[10 350 120 20]);
agentsBox = uieditfield(ctrl2,'numeric','Value',50,'Position',[140 350 80 20]);
uilabel(ctrl2,'Text','Number of Trials:','Position',[10 310 120 20]);
trialBox2 = uieditfield(ctrl2,'numeric','Value',300,'Position',[140 310 80 20]);

uibutton(ctrl2,'Text','Run Simulation','Position',[50 260 150 40],...
    'ButtonPushedFcn',@runMultiAgent);

% Axes
axQ2 = uiaxes(tab2,'Position',[300 300 680 200]);
title(axQ2,'Average Q-values'); xlabel(axQ2,'Trial'); ylabel(axQ2,'Q'); hold(axQ2,'on'); grid(axQ2,'on');
axF2 = uiaxes(tab2,'Position',[300 50 680 200]);
title(axF2,'Strategy Frequencies'); xlabel(axF2,'Trial'); ylabel(axF2,'Fraction'); hold(axF2,'on'); grid(axF2,'on');

%% ---------------- TAB 3: BEST ROLE TRACKING -----------------
tab3 = uitab(tg,'Title','Best Role Evaluation');
uilabel(tab3,'Text',['Multiple learners explore all roles to estimate which role is best. ',...
    'Shows cumulative rewards per role over time.'],...
    'Position',[20 500 950 50],'FontSize',14);

ctrl3 = uipanel(tab3,'Title','Controls','Position',[10 10 250 480]);
uilabel(ctrl3,'Text','Learning Rate α:','Position',[10 430 120 20]);
lrBox3 = uieditfield(ctrl3,'numeric','Value',0.2,'Limits',[0 1],'Position',[140 430 80 20]);
uilabel(ctrl3,'Text','Softmax τ:','Position',[10 390 120 20]);
tauBox3 = uieditfield(ctrl3,'numeric','Value',1,'Position',[140 390 80 20]);
uilabel(ctrl3,'Text','Number of Learners:','Position',[10 350 120 20]);
learnersBox = uieditfield(ctrl3,'numeric','Value',50,'Position',[140 350 80 20]);
uilabel(ctrl3,'Text','Number of Trials:','Position',[10 310 120 20]);
trialBox3 = uieditfield(ctrl3,'numeric','Value',300,'Position',[140 310 80 20]);

uibutton(ctrl3,'Text','Run Simulation','Position',[50 260 150 40],...
    'ButtonPushedFcn',@runBestRole);

% Axes
axCum3 = uiaxes(tab3,'Position',[300 50 680 450]);
title(axCum3,'Cumulative Reward per Role'); xlabel(axCum3,'Trial'); ylabel(axCum3,'Cumulative Reward'); hold(axCum3,'on'); grid(axCum3,'on');

%% ---------------- COLORS & STRATEGY NAMES -----------------
strategies = {'Hawk','Dove','Retaliator','Assessor'};
colors = [0 0 0.7; 0 0.6 0; 0.85 0.33 0; 0.8 0.6 0.1]; % Colorblind friendly

%% ---------------- SIMULATION FUNCTIONS -----------------
    function runSingleLearner(~,~)
        alpha = lrBox1.Value; tau = tauBox1.Value; T = trialBox1.Value;
        Q = zeros(1,4); Q_hist = zeros(T,4); P_hist = zeros(T,4);
        V = 50; C = 100;
        M = [ (V-C)/2,  V,      (V-C)/2, V;
              0,       V/2,     V/2,    0;
              (V-C)/2, V/2,     V/2,    V/2;
              (V-C)/2, V,       V/2,    V/2];
        for t = 1:T
            expQ = exp(Q/tau); policy = expQ/sum(expQ);
            a = randsample(1:4,1,true,policy); opp = randi(4);
            r = M(a,opp);
            Q(a) = Q(a)+alpha*(r-Q(a));
            Q_hist(t,:) = Q; P_hist(t,:) = policy;
        end
        cla(axQ1); cla(axP1);
        for i=1:4
            plot(axQ1,Q_hist(:,i),'Color',colors(i,:),'LineWidth',2); hold(axQ1,'on');
            plot(axP1,P_hist(:,i),'Color',colors(i,:),'LineWidth',2); hold(axP1,'on');
        end
        legend(axQ1,strategies,'Location','northwest');
        legend(axP1,strategies,'Location','northwest');
    end

    function runMultiAgent(~,~)
        alpha = lrBox2.Value; tau = tauBox2.Value;
        N = agentsBox.Value; T = trialBox2.Value;
        V = 50; C = 100;
        M = [ (V-C)/2,  V,      (V-C)/2, V;
              0,       V/2,     V/2,    0;
              (V-C)/2, V/2,     V/2,    V/2;
              (V-C)/2, V,       V/2,    V/2];
        Q = zeros(N,4); % each agent
        Q_hist = zeros(T,4); F_hist = zeros(T,4);
        for t = 1:T
            actions = zeros(N,1);
            rewards = zeros(N,1);
            for i = 1:N
                expQ = exp(Q(i,:)/tau); policy = expQ/sum(expQ);
                a = randsample(1:4,1,true,policy); opp = randsample(1:4,1);
                r = M(a,opp); Q(i,a) = Q(i,a)+alpha*(r-Q(i,a));
                actions(i) = a; rewards(i) = r;
            end
            Q_hist(t,:) = mean(Q,1);
            for k = 1:4
                F_hist(t,k) = sum(actions==k)/N;
            end
        end
        cla(axQ2); cla(axF2);
        for i=1:4
            plot(axQ2,Q_hist(:,i),'Color',colors(i,:),'LineWidth',2); hold(axQ2,'on');
            plot(axF2,F_hist(:,i),'Color',colors(i,:),'LineWidth',2); hold(axF2,'on');
        end
        legend(axQ2,strategies,'Location','northwest');
        legend(axF2,strategies,'Location','northwest');
    end

    function runBestRole(~,~)
        alpha = lrBox3.Value; tau = tauBox3.Value; N = learnersBox.Value; T = trialBox3.Value;
        V = 50; C = 100;
        M = [ (V-C)/2,  V,      (V-C)/2, V;
              0,       V/2,     V/2,    0;
              (V-C)/2, V/2,     V/2,    V/2;
              (V-C)/2, V,       V/2,    V/2];
        Q = zeros(N,4); cumRewards = zeros(1,4); cumHist = zeros(T,4);
        for t = 1:T
            actions = zeros(N,1); rewards = zeros(N,1);
            for i = 1:N
                expQ = exp(Q(i,:)/tau); policy = expQ/sum(expQ);
                a = randsample(1:4,1,true,policy); opp = randsample(1:4,1);
                r = M(a,opp); Q(i,a)=Q(i,a)+alpha*(r-Q(i,a));
                cumRewards(a)=cumRewards(a)+r; actions(i)=a;
            end
            cumHist(t,:) = cumRewards;
        end
        cla(axCum3);
        for i=1:4
            plot(axCum3,cumHist(:,i),'Color',colors(i,:),'LineWidth',2); hold(axCum3,'on');
        end
        legend(axCum3,strategies,'Location','northwest');
    end

end



