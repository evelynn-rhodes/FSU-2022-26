function indiv_simulation_window
% INDIV_SIMULATION_WINDOW
% Unified GUI + engine for INDIVIDUAL Hawk–Dove–Retaliator–Assessor simulation

    f = uifigure('Name','Individual Dynamics Simulation','Position',[200 200 650 500]);

    % -----------------------------------------------------
    % Tab group
    % -----------------------------------------------------
    t = uitabgroup(f,'Position',[10 10 630 480]);

    tab1 = uitab(t,'Title','Variable Description');
    tab2 = uitab(t,'Title','Simulation');

    %% --------------------- TAB 1: EXPLANATION -----------------------------
    uitextarea(tab1, ...
        'Parent', tab1, ...
        'Position',[10 10 600 440], ...
        'Editable','off', ...
        'Value', [
            "INDIVIDUAL DYNAMICS"
            "----------------------"
           "This individual dynamics simulation models a single encounter between two pigeons, with payoffs determined by a 4×4 matrix M that quantifies the net fitness consequences of each possible role interaction. Here, V is the value of winning, set in other simulations to 50, and C is the cost of fighting, set to 100; these two numbers define how much net benefit or loss each interaction produces. The payoff matrix M encodes all outcomes: the rows correspond to the strategy of Player 1 (Hawk, Dove, Retaliator, Assessor) and the columns correspond to Player 2’s strategy. Each entry M(i,j) represents the expected net fitness Player 1 receives when using strategy i against strategy j, accounting for both gains from winning and losses from conflict. The simulation allows the user to enter a custom M, so they can adjust how severe the costs of fights are relative to the rewards, and then select two roles to interact. When the simulation runs, it looks up the appropriate payoffs for each player in M, narrates the behavior and outcome, and shows who gains the advantage. This setup emphasizes the logic of game-theoretic interactions: V and C determine the relative incentives for aggression versus peace, and M translates those incentives into explicit fitness units that drive decision-making in a single encounter."
        ]);

    %% --------------------- TAB 2: SIMULATION -----------------------------
    % Labels
    uilabel(tab2,'Parent',tab2,'Text','Player 1 Strategy:','Position',[40 420 150 22]);
    dd1 = uidropdown(tab2,'Parent',tab2,'Items',{'Hawk','Dove','Retaliator','Assessor'}, ...
                     'Position',[180 420 130 22]);

    uilabel(tab2,'Parent',tab2,'Text','Player 2 Strategy:','Position',[40 380 150 22]);
    dd2 = uidropdown(tab2,'Parent',tab2,'Items',{'Hawk','Dove','Retaliator','Assessor'}, ...
                     'Position',[180 380 130 22]);

    uilabel(tab2,'Parent',tab2,'Text','4×4 Payoff Matrix:','Position',[40 320 150 22]);

    matrixBox = uitextarea(tab2, ...
        'Parent',tab2, ...
        'Value', {
            '0   6   2   3'
            '-2  3   3   2'
            '1   2   4   3'
            '1   1   2   5'
        }, ...
        'Position',[40 200 300 120]);

    % Run button
    runBtn = uibutton(tab2,'Parent',tab2,'Text','Run Simulation', ...
        'Position',[380 380 150 40], ...
        'ButtonPushedFcn', @(~,~) runSim() );

    % Output area
    narrBox = uitextarea(tab2, ...
        'Parent',tab2, ...
        'Position',[40 20 580 160], ...
        'Editable','off');

    %% --------------------- CALLBACK -------------------------------------
    function runSim()
        try
            M = str2num(char(matrixBox.Value)); %#ok<ST2NM>
            if ~isequal(size(M),[4 4])
                narrBox.Value = "Error: Payoff matrix must be 4×4.";
                return;
            end

            s1 = dd1.Value;
            s2 = dd2.Value;

            [narr, ~] = hawkdove_engine(s1, s2, M);

            narrBox.Value = narr;
            drawnow;

        catch ME
            narrBox.Value = sprintf("Error:\n%s", ME.message);
        end
    end

    %% =====================================================================
    %%                 NESTED ENGINE FUNCTION
    %% =====================================================================
    function [narration, payoff] = hawkdove_engine(s1, s2, M)

        roles = {'Hawk','Dove','Retaliator','Assessor'};
        idx = @(s) find(strcmpi(s, roles));

        i = idx(s1);
        j = idx(s2);

        p1 = M(i,j);
        p2 = M(j,i);
        payoff = [p1 p2];

        % ---- Behavior descriptions ---------------------------------------
        switch s1
            case 'Hawk'
                act1 = "Player 1 acts aggressively and commits to fighting.";
            case 'Dove'
                act1 = "Player 1 avoids escalation and attempts peaceful display.";
            case 'Retaliator'
                act1 = "Player 1 stays peaceful unless attacked, then retaliates.";
            case 'Assessor'
                act1 = "Player 1 evaluates opponent strength before deciding to escalate.";
        end

        switch s2
            case 'Hawk'
                act2 = "Player 2 escalates aggressively into a fight.";
            case 'Dove'
                act2 = "Player 2 displays peacefully without escalation.";
            case 'Retaliator'
                act2 = "Player 2 is peaceful unless attacked, then fights back.";
            case 'Assessor'
                act2 = "Player 2 assesses whether to escalate based on strength.";
        end

        % ---- Winner/loser narration --------------------------------------
        if p1 > p2
            resultTxt = "Player 1 gains the advantage and receives the higher payoff.";
        elseif p2 > p1
            resultTxt = "Player 2 gains the advantage and receives the higher payoff.";
        else
            resultTxt = "Both players receive the same payoff from this encounter.";
        end

        narration = sprintf( ...
            "Interaction Summary:\n" + ...
            "\nPlayer 1 (%s):\n  %s" + ...
            "\n\nPlayer 2 (%s):\n  %s" + ...
            "\n\nOutcome:\n  Player 1 payoff: %g units\n  Player 2 payoff: %g units\n\n%s", ...
            s1, act1, ...
            s2, act2, ...
            p1, p2, resultTxt );

    end % end nested engine

end
