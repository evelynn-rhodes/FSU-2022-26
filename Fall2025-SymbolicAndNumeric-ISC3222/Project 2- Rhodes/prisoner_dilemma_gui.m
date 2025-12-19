function prisoner_dilemma_gui
    % PRISONER_DILEMMA_GUI
    % Fully functioning educational Prisoner's Dilemma simulation
    % Uses jail years instead of payoff points.

    %% Create main figure
    fig = uifigure('Name','Prisoner''s Dilemma','Position',[200 200 800 500]);

    %% Tab group
    tg = uitabgroup(fig,'Position',[0 0 fig.Position(3) fig.Position(4)]);

    % ---------------- Single Player Tab ----------------
    tab1 = uitab(tg,'Title','Single Player');

    % Instructions
    uitextarea(tab1,'Value',[
        "Welcome to the Prisoner's Dilemma!"
        "---------------------------------"
        "You and another gang member have been arrested and separated. Each of you must choose between staying silent (Cooperate) or testifying for a lighter sentence (Defect)."
        ""
        "Jail Time Rules:"
        "• If both stay silent: 1 year each"
        "• If you stay silent and they testify: you get 3 years, they get 0"
        "• If you testify and they stay silent: you get 0, they get 3"
        "• If both testify: 5 years each"
        ""
        "Choose your move each round."
        "The computer can use several strategies. Altruistic means the comptuer will always stay silent, defector means the comptuer will always testify, and tit-for-tat is the computer will testify next round if you testified the last round. Random is a coin flip of if the comptuer will testify or stay silent."
    ],'Editable',false,'Position',[20 250 350 200]);

    % User buttons
    btnCoop = uibutton(tab1,'Text','Cooperate','Position',[400 350 100 40]);
    btnDef  = uibutton(tab1,'Text','Defect','Position',[520 350 100 40]);

    % Computer strategy dropdown
    uilabel(tab1,'Text','Computer Strategy','Position',[400 300 120 20]);
    ddStrat = uidropdown(tab1,'Items',{'Altruistic','Defector','Tit-for-Tat','Random'},...
        'Position',[520 300 100 25]);

    % Outcome axes
    axOutcome = axes(tab1,'Position',[0.5 0.05 0.45 0.4]);
    axis(axOutcome,[0 2 0 1]);
    axOutcome.XTick = []; axOutcome.YTick = []; 
    hold(axOutcome,'on');
    title(axOutcome,'Outcome: You vs Opponent');

    % Output text
    txtDesc = uitextarea(tab1,'Position',[20 50 350 180],'Editable',false,'Value',{''});

    % Variables
    userScore = 0;
    compScore = 0;
    lastUserMove = 1;
    roundNum = 0;

    % Callbacks
    btnCoop.ButtonPushedFcn = @(~,~)playRound(1);
    btnDef.ButtonPushedFcn  = @(~,~)playRound(2);

    function playRound(userMove)
        roundNum = roundNum + 1;
        compStrategy = ddStrat.Value;

        % Determine computer move
        switch compStrategy
            case 'Altruistic'
                compMove = 1;
            case 'Defector'
                compMove = 2;
            case 'Tit-for-Tat'
                if roundNum==1
                    compMove = 1;
                else
                    compMove = lastUserMove;
                end
            case 'Random'
                compMove = randi(2);
        end

        % Jail time matrix (years)
        % Rows: user [C=1, D=2]
        % Cols: comp [C=1, D=2]
        J = [1 3;   % user cooperates
             0 5];  % user defects

        % Assign jail years
        userYears = J(userMove, compMove);
        compYears = J(compMove, userMove);

        % Update cumulative totals
        userScore = userScore + userYears;
        compScore = compScore + compYears;

        % Outcome visuals
        cla(axOutcome)
        if userMove==1 && compMove==1
            rectangle(axOutcome,'Position',[0.2 0.2 0.3 0.6],'FaceColor','g')
            rectangle(axOutcome,'Position',[1   0.2 0.3 0.6],'FaceColor','g')
        elseif userMove==1 && compMove==2
            rectangle(axOutcome,'Position',[0.2 0.2 0.3 0.6],'FaceColor','r')
            rectangle(axOutcome,'Position',[1   0.2 0.3 0.6],'FaceColor','g')
        elseif userMove==2 && compMove==1
            rectangle(axOutcome,'Position',[0.2 0.2 0.3 0.6],'FaceColor','g')
            rectangle(axOutcome,'Position',[1   0.2 0.3 0.6],'FaceColor','r')
        else
            rectangle(axOutcome,'Position',[0.2 0.2 0.3 0.6],'FaceColor','r')
            rectangle(axOutcome,'Position',[1   0.2 0.3 0.6],'FaceColor','r')
        end

        % Update text description
        moves = {'Cooperate','Defect'};
        txtDesc.Value = {...
            sprintf('Round %d',roundNum),...
            sprintf('You chose: %s',moves{userMove}),...
            sprintf('Opponent chose: %s',moves{compMove}),...
            sprintf('Years in jail this round: You = %d, Opponent = %d',userYears,compYears),...
            sprintf('Total years in jail: You = %d, Opponent = %d',userScore,compScore)
        };

        lastUserMove = userMove;
    end

    % ---------------- Two-Player Tab ----------------
    tab2 = uitab(tg,'Title','Two Player');

    uitextarea(tab2,'Value',[
        "Two-Player Prisoner's Dilemma!"
        "-------------------------------"
        "You and another gang member have been arrested and separated. Each of you must choose between staying silent (Cooperate) or testifying for a lighter sentence (Defect)."
        "Player 1 and Player 2 choose in secret. After both moves, jail years are shown."
        ""
        "Jail Time Rules:"
        "• If both stay silent: 1 year each"
        "• If you stay silent and they testify: you get 3 years, they get 0"
        "• If you testify and they stay silent: you get 0, they get 3"
        "• If both testify: 5 years each"
        ""
    ],'Editable',false,'Position',[20 250 350 200]);

    % Buttons
    btnP1C = uibutton(tab2,'Text','P1 Cooperate','Position',[400 350 100 40]);
    btnP1D = uibutton(tab2,'Text','P1 Defect','Position',[520 350 100 40]);
    btnP2C = uibutton(tab2,'Text','P2 Cooperate','Position',[400 250 100 40],'Enable','off');
    btnP2D = uibutton(tab2,'Text','P2 Defect','Position',[520 250 100 40],'Enable','off');

    % Outcome axes
    axOutcome2 = axes(tab2,'Position',[0.5 0.05 0.45 0.4]);
    axis(axOutcome2,[0 2 0 1]); 
    hold(axOutcome2,'on'); 
    axOutcome2.XTick = []; axOutcome2.YTick = [];
    title(axOutcome2,'Outcome: Player 1 vs Player 2');

    txtDesc2 = uitextarea(tab2,'Position',[20 50 350 180],'Editable',false,'Value',{''});

    % Stored moves
    movesP = [0 0];

    % Callbacks
    btnP1C.ButtonPushedFcn = @(~,~)selectMove(1,1);
    btnP1D.ButtonPushedFcn = @(~,~)selectMove(1,2);
    btnP2C.ButtonPushedFcn = @(~,~)selectMove(2,1);
    btnP2D.ButtonPushedFcn = @(~,~)selectMove(2,2);

    function selectMove(player,choice)
        movesP(player) = choice;

        if player==1
            btnP1C.Enable = 'off'; btnP1D.Enable = 'off';
            btnP2C.Enable = 'on';  btnP2D.Enable = 'on';
        else
            evaluateTwoPlayerRound();
            btnP1C.Enable = 'on'; btnP1D.Enable = 'on';
            btnP2C.Enable = 'off'; btnP2D.Enable = 'off';
            movesP = [0 0];
        end
    end

    function evaluateTwoPlayerRound()
        % Jail-time matrix
        J = [1 3;
             0 5];

        P1Years = J(movesP(1), movesP(2));
        P2Years = J(movesP(2), movesP(1));

        % Graphics
        cla(axOutcome2)
        if movesP(1)==1 && movesP(2)==1
            rectangle(axOutcome2,'Position',[0.2 0.2 0.3 0.6],'FaceColor','g')
            rectangle(axOutcome2,'Position',[1   0.2 0.3 0.6],'FaceColor','g')
        elseif movesP(1)==1 && movesP(2)==2
            rectangle(axOutcome2,'Position',[0.2 0.2 0.3 0.6],'FaceColor','r')
            rectangle(axOutcome2,'Position',[1   0.2 0.3 0.6],'FaceColor','g')
        elseif movesP(1)==2 && movesP(2)==1
            rectangle(axOutcome2,'Position',[0.2 0.2 0.3 0.6],'FaceColor','g')
            rectangle(axOutcome2,'Position',[1   0.2 0.3 0.6],'FaceColor','r')
        else
            rectangle(axOutcome2,'Position',[0.2 0.2 0.3 0.6],'FaceColor','r')
            rectangle(axOutcome2,'Position',[1   0.2 0.3 0.6],'FaceColor','r')
        end

        % Description
        moves = {'Cooperate','Defect'};
        txtDesc2.Value = {...
            sprintf('Player 1 chose: %s',moves{movesP(1)}),...
            sprintf('Player 2 chose: %s',moves{movesP(2)}),...
            sprintf('Years in jail: P1 = %d, P2 = %d',P1Years,P2Years)
        };
    end
end

