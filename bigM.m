classdef bigM < handle
    % BIGM This class implements big M Simplex Method to solve a
    % linear programming problem in the following format
    % min/max c'x
    % s.t.   Ax {>=, =, <=} b
    % x >= 0
    %
    % The code here follows the variation taught in the FuO Course for 
    % OR Msc. Specifically, this class is designed for class demonstration 
    % and small problems only, neither suitable for large problem nor 
    % powerful enough for high performance computing.
    %
    % Example: See example.m
    %
    % 3rd Oct 2013
    % Yiming YAN 
    % University of Edinburgh
    
    properties(SetAccess = private)
        A;                      % Matrix in Ax {<=, =, >=} b
        b;                      % RHS
        c;                      % Vector of coefficients of objective function 
        inq;                    % m dimentional vector, 
                                % 1 for >=; 0 for ==; -1 for <=.
        m;                      % Number of constraints
        n;                      % Number of original primal variables
        x_opt;                  % Optimal solution
        fval;                   % Optimal objective value
        type;                   % Minimise or maximise the problem
    end
    
    properties( Access = private )
        M;                      % bigM
        
        n_slack;                % Number of slack vars added
        n_artfVar;              % number of artificial vars added
        artfcl_var_set          % The index set of artificial variables
        artfcl_var_removed = 1; % 1, if removed
        
        counter = 1;            % Iteration counter
        
        basis;                  % Vector of indices of basic basic variables
        reducedCost;             
        whichRow;                
        
        tab;                    % Tableau
        status;                 % maxIter, unbounded, infeasible, optimal
        terminate;
        vars_list;              % List of names of all variables
        unbounded_vars;         % unbounded variables
    end
    
    properties( Constant )
        maxIter = 20;
    end
    
    methods
        
        function p = bigM(A,b,c,inq, type)
            % CONSTRUCTOR
            
            % Check input
            if nargin < 4
                error('bigM: Not enough input');
            end
            
            %% Read the input
            % Check min or max
            if strcmpi(type,'max')
                c = -c;
            end
            
            % Check if b is negative
            b_Neg = b < 0;
            inq(b_Neg) = -inq(b_Neg);
            b(b_Neg) = -b(b_Neg);
            A(b_Neg,:) = -A(b_Neg,:);
            
            % Assiqn properties
            p.A = A; p.b = b; p.c = c;
            p.inq = inq;
            p.terminate = 0;
            
            % Set the value of big M!
            p.M = ceil(max(norm(A),norm(b))/100)*100;
            
            [p.m, p.n] = size(A);
            
            % Form the vars_list
            p.vars_list = cell(1,p.n);
            for i=1:p.n
                p.vars_list{i} = ['x' num2str(i)];
            end
            
            p.type = type;
        end
        
        %% Driver function
        function solve(p)
            % SOLVE This is the driver function for bigM which actually
            % solves the given problem.
            
            p.transform_to_standardForm;
            p.construct_initial_tableau;
            p.printTab;
            while ~p.terminate
                p.do_simplex_iterate;
                p.increment_counter;
                p.printTab;
                p.remove_artfcl_vars;
                p.check_termination;
            end
            p.output_summary;
        end
        
        function transform_to_standardForm(p)
            
            % add slack variables
            p.n_slack = sum(p.inq < 0)+sum(p.inq > 0);
            p.A = [p.A zeros(p.m, p.n_slack)];
            p.c = [p.c; zeros(p.n_slack,1)];
            
            idx_slack = p.n+1;
            for i = 1:p.m
                switch p.inq(i)
                    case 1  % >=
                        p.A(i, idx_slack) = -1;
                        idx_slack = idx_slack+1;
                         p.vars_list{end+1} = ['s' num2str(i)];
                    case -1 % <=
                        p.A(i, idx_slack) = 1;
                        idx_slack = idx_slack+1;
                        p.vars_list{end+1} = ['s' num2str(i)];
                end
            end
        end
        
        function construct_initial_tableau(p)

            % Find potential basis
            for i = 1:p.n+p.n_slack
                
                if nnz(p.A(:,i)) == 1
                    row_number = find(p.A(:,i) == 1);
                    if ~isempty(row_number)                        
                        % Add the col if the current row has not been selected
                        if ~ismember(row_number, p.whichRow)
                            p.whichRow(end+1) = row_number;
                            p.basis(row_number) = i;
                        end
                    end
                end
            end
            
            % Add bigM if necessary            
            p.n_artfVar = p.m - sum(p.basis > 0);
            if  p.n_artfVar > 0
                p.artfcl_var_removed = 0;
                
                % Record the index for artificial variables
                p.artfcl_var_set = (p.n + p.n_slack + 1) : (p.n + p.n_slack + p.n_artfVar);
                
                p.A = [p.A zeros(p.m,p.n_artfVar)];
                
                add_to_rows = setdiff(1:p.m, p.whichRow);
                
                % Formulate matrix A with newly added artificial vars
                for i = 1:length(add_to_rows)
                   p.A(add_to_rows(i), p.n+p.n_slack+i) = 1;
                   p.vars_list{end+1} = ['a' num2str(i)];
                   p.basis(add_to_rows(i)) = p.n+p.n_slack+i;
                end
                
                p.c = [p.c; p.M*ones(p.n_artfVar,1)];
                
            end
                        
            % Get reducedCost
            p.reducedCost = zeros(1,p.n + p.n_slack + p.n_artfVar);
            y = p.c(p.basis);
            
            p.reducedCost = (p.c - p.A'*y)';
            
            % Construct initial tableau
            p.tab = [p.A  p.b];
            p.tab = [p.tab; [ p.reducedCost -p.b'*y ]];
        end
        
        function do_simplex_iterate(p)
            % DO_SIMPLEX_ITERATE Perform one iterate of simplex method.
            
            % Find the smallest reduced cost
            [~, workCol] = min( p.tab(end,1:end-1) );
            
            % Choose the smallest postive ratio
            work_col_positive_elems = find( p.tab(1:end-1,workCol) > 0 );
            
            [~, chooseRow] = min( p.tab(work_col_positive_elems, end) ./...
                p.tab(work_col_positive_elems, workCol) );
            chooseRow = work_col_positive_elems(chooseRow);
            
            % In/out basis
            p.basis(chooseRow) = workCol;
            
            % Get pivot
            pivot = p.tab(chooseRow, workCol);
            
            % Elemenate elements in work cols
            % and update the value of elements in the tableau.
            p.tab(chooseRow, :) = p.tab(chooseRow,:) / pivot;
            
            tmp_row_indx = setdiff(1:p.m+1,chooseRow);
            p.tab(tmp_row_indx, :) = p.tab(tmp_row_indx, :) -...
                ones(p.m,1) * p.tab(chooseRow,:) .*...
                ( p.tab(tmp_row_indx, workCol) *...
                ones(1,size(p.tab,2)) );
        end
        
        function remove_artfcl_vars(p)
            % Check if artificial variables can be removed
            if ~p.artfcl_var_removed
                if all( p.tab(end, p.artfcl_var_set) > 0 )
                    p.tab(:,p.artfcl_var_set) = [];
                    p.vars_list( p.artfcl_var_set ) = [];
                    p.artfcl_var_removed = 1;
                    p.printTab;
                end
            end
        end
        
        function check_termination(p)
            p.terminate = 0;
            
            % Reach maxIter
            if p.counter >= p.maxIter
                p.terminate = 1;
                p.status = 'maxIter';
            end
            
            % Nonnegative reduced costs - optimal
            if all(p.tab(end,1:end-1) > -1e-10)
                p.terminate = 1;
                
                % infeasible
                if ~p.artfcl_var_removed
                    p.status = 'infeasible';
                else
                    % optimal
                    p.status = 'optimal';
                    p.x_opt = zeros(p.n + p.n_slack + p.n_artfVar, 1);
                    p.x_opt(p.basis) = p.tab(1:end-1,end);
                    if strcmpi(p.type,'max')
                        p.fval = p.tab(end,end);
                    elseif strcmpi(p.type,'min')
                        p.fval = -p.tab(end,end);
                    end
                end
            else
                % unbounded
                negative_rc_cols = find(p.tab (end,1:end-1) < 0);
                unbounded_idx = sum(p.tab(1:end-1, negative_rc_cols) < 0);
                unbounded_idx = unbounded_idx == p.m;
                if any(unbounded_idx)
                    p.status = 'unbounded';
                    p.terminate = 1;
                    p.unbounded_vars = negative_rc_cols(unbounded_idx > 0);
                end                
            end
        end
        
        
        function increment_counter(p)
            p.counter = p.counter+1;
        end
        
        function printTab(p)
            fprintf('\n\n');
            fprintf('\t\tTableau %2s:\n\n', num2str(p.counter));
            % Print header
            for k=1:length(p.vars_list)+1
                if k==1
                    fprintf('%2s | ', '');
                else
                    fprintf('%8s %4s', p.vars_list{k-1}, '');
                    if k == length(p.vars_list)+1
                        fprintf(' | %8s', 'RHS');
                    end
                end
            end
            fprintf('\n');
            fprintf('------------------------------------------------------------------\n');
            
            for i=1:p.m+1
                if i <= p.m
                    fprintf('%2s | ', p.vars_list{p.basis(i)});
                else
                    fprintf('------------------------------------------------------------------\n');
                    fprintf('%2s | ', 'RC');
                end
                for j = 1:size(p.tab,2)
                    fprintf('%8.2f %4s', p.tab(i,j), '');
                    if j == size(p.tab,2)-1
                       fprintf(' | ') 
                    end
                end
                fprintf('\n');
            end
            
        end
        
        function output_summary(p)
            if p.terminate
                fprintf('\n========== ========== ==========\n');
                fprintf('Terminated. \n');
                switch lower( p.status )
                    case 'optimal'
                        fprintf('Optimal solution: \n');
                        for i=1:p.n
                            fprintf('x[%2d ] = %5.2f\n', i, p.x_opt(i));
                        end
                        fprintf('Optimal objective function value: %5.2f\n', p.fval);
                        
                    case 'unbounded'
                        fprintf('This problem is unbounded. \n');
                        fprintf('Unbounded variables: ');
                        
                        for i = 1:length(p.unbounded_vars)
                            fprintf('%2s ',p.vars_list{p.unbounded_vars(i)});
                        end
                        fprintf('\n');
                        
                    case 'infeasible'
                        fprintf('This problem is infeasible.\n');
                    case 'maxiter'
                        fprintf('Maximum number of iterations reached.\n');
                end
                
            end
        end
    end
end