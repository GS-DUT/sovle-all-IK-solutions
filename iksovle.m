function [qt,tcount,e0] = iksovle(robot, tr, q, varargin)
    
    n = robot.n;
    
    TT = SE3.check(tr);
    
    %  set default parameters for solution
    opt.ilimit = 50;
    opt.rlimit = 20;
    opt.slimit = 100;
    opt.tol = 1e-10;
    opt.lambda = 0.0001;
    opt.lambdamin = 1e-6;
    opt.search = false;
    opt.quiet = false;
    opt.verbose = false;
    opt.mask = [1 1 1 1 1 1];
    opt.q0 = zeros(1, n);
    opt.transpose = NaN;
    
    [opt,args] = tb_optparse(opt, varargin);
    
%     qlim = zeros(n,2);
    
    assert(numel(opt.mask) == 6, 'RTB:ikine:badarg', 'Mask matrix should have 6 elements');
    assert(n >= numel(find(opt.mask)), 'RTB:ikine:badarg', 'Number of robot DOF must be >= the same number of 1s in the mask matrix');
    W = diag(opt.mask);
    
    
    qt = zeros(length(TT), n);  % preallocate space for results 
    tcount = 0;              % total iteration count
    
    failed = false;
    revolutes = robot.isrevolute();
    
    e0 = [];
    
    for i=1:length(TT)
        T = TT(i);
        lambda = opt.lambda;

        if i>1
            q = qt(i-1,:);
        end

        iterations = 0;
        
        if opt.debug
            e = tr2delta(robot.fkine(q), T);
            fprintf('Initial:  |e|=%g\n', norm(W*e));
        end
        
        while true
            % update the count and test against iteration limit
            iterations = iterations + 1;
            if iterations > opt.ilimit
                if ~opt.quiet
                    warning('ikine: iteration limit %d exceeded (pose %d), final err %g', ...
                        opt.ilimit, i, nm);
                end
                failed = true;
                break
            end
            
            e = tr2delta(robot.fkine(q), T);
            e0(iterations) = norm(e);
            % are we there yet
            if norm(W*e) < opt.tol  
                break;
            end
            
            % compute the Jacobian
            J = jacobe(robot, q+[0,pi/2,0,0,0,pi/2]);
            
            JtJ = J'*W*J;
            
            u = 0.001*norm(e);
            
            if ~isnan(opt.transpose)
                % do the simple Jacobian transpose with constant gain
                dq = opt.transpose * J' * e;
            else
                % do the damped inverse Gauss-Newton with Levenberg-Marquadt
                dq = inv(JtJ + u * eye(size(JtJ)) ) * J' * W * e;
                % compute possible new value of
                qnew = q + dq';

                % and figure out the new error
                enew = tr2delta(robot.fkine(qnew), T);

%                 r = -(norm(W*e)^2-norm(W*enew)^2)/(norm(W*e)^2-norm(W*(e+J*dq))^2);
                
                    if opt.debug
                        fprintf('ACCEPTED: |e|=%g, |dq|=%g, lambda=%g\n', norm(W*enew), norm(dq), lambda);
                    end
                    q = qnew;
                    e = enew;
                    
                    % step is rejected, increase the damping and retry
                    if opt.debug
                        fprintf('rejected: |e|=%g, |dq|=%g, lambda=%g\n', norm(W*enew), norm(dq), lambda);
                    end
                    
            end
            
            
            % wrap angles for revolute joints
            k = (q > pi) & revolutes;
            q(k) = q(k) - 2*pi;
            
            k = (q < -pi) & revolutes;
            q(k) = q(k) + 2*pi;
            
            nm = norm(W*e);
            
            
        end  % end ikine solution for this pose
        qt(i,:) = q';
        tcount = tcount + iterations;
        if opt.verbose && ~failed
            fprintf('%d iterations\n', iterations);
        end
        if failed
            if ~opt.quiet
                warning('failed to converge: try a different initial value of joint coordinates');
            end
%             qt = [];
            qt = zeros(length(TT), n);
        end
    end
    
    
    if opt.verbose && length(TT) > 1
        fprintf('TOTAL %d iterations\n', tcount);
    end
    
    
end
