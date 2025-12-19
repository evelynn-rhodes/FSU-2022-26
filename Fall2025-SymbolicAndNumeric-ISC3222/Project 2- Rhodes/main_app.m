function main_app
% MAIN_APP
% Top-level GUI for the educational Game Theory app

    f = uifigure('Name','Game Theory Explorer','Position',[200 200 600 500]);

    %% OVERVIEW TEXT
    overview = [
        "GAME THEORY OVERVIEW"
        "---------------------------"
        ""
        "Game theory is the study of decision-making in situations where the outcome for each participant depends not only on their own choices but also on the choices of others. At its core, game theory is about strategies, payoffs, and interactions. A strategy is simply a plan or action that a person (or “player”) can take. The payoff is the reward or consequence associated with a combination of strategies, which can represent money, time, happiness, fitness, or any measurable outcome. The simplest building blocks are players, actions, and payoffs, and from these, more complex ideas emerge, like Nash equilibria—situations where no player can benefit by changing their strategy alone—and cooperation versus competition. Game theory applies everywhere in everyday life: deciding whether to share responsibilities at work, negotiating prices, competing for limited resources, or even choosing whether to wait your turn in line. It’s also foundational in biology for understanding animal behavior, in economics for markets, in politics for strategy, and in computer science for algorithms. By modeling choices and consequences systematically, game theory gives a framework for predicting outcomes and understanding the incentives that drive behavior."
        ""
        "Click below to explore the examples of game theory."
    ];

    uitextarea(f,'Value',overview,'Editable','off',...
        'Position',[20 100 560 360]);

    %% BUTTON TO ENTER PIGEON MODULE
    uibutton(f,'Text','Pigeon Example','Position',[20 40 200 40], ...
        'ButtonPushedFcn', @(~,~) openPigeonMenu());
    uibutton(f,'Text','Prisoners dilema','Position',[350 40 200 40], ...
        'ButtonPushedFcn', @(~,~) prisoner_dilemma_gui());

    %% ------------------ NESTED FUNCTION -----------------------
    function openPigeonMenu()
        fig = uifigure('Name','Pigeon Game Theory','Position',[220 220 620 520]);

        txt = {
            'THE PIGEON (HAWK–DOVE) MODEL'
            '--------------------------------------------------------------'
            ''
            'This module describes the classic Maynard Smith model of animal'
            'conflict. Each strategy represents a rule for escalation or retreat:'
            ''
            '  • Hawk       — escalates into a fight'
            '  • Dove       — avoids escalation and shares or retreats'
            '  • Retaliator — Dove unless attacked, then fights'
            '  • Assessor   — fights only when stronger'
            ''
            'Choose a simulation:'
        };

        uitextarea(fig,'Value',txt,'Editable','off',...
            'Position',[20 180 580 300]);

        uibutton(fig,'Text','Individual Dynamics', ...
            'Position',[50 80 150 40], ...
            'ButtonPushedFcn', @(~,~) indiv_simulation_window());
        uibutton(fig,'Text','Learning Simulation', ...
            'Position',[250 80 150 40], ...
            'ButtonPushedFcn', @(~,~) learning_simulation_window());
        uibutton(fig,'Text','Population Dynamics', ...
            'Position',[450 80 150 40], ...
            'ButtonPushedFcn', @(~,~) pop_dynamics_window());
    end

end

