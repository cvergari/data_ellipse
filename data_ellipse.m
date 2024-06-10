classdef data_ellipse < handle
    % data_ellipse is an object that allows to calculate a data ellipse,
    % containing 95% of two linearly related and normally distributed
    % varaibles.
    %
    % To see it in action, type :
    %   obj = data_ellipse;
    %   obj.demo()
    %   
    % Example:
    %   obj = data_ellipse;
    %
    %   % Make linearly correlated synthetic data 
    %   obj.x_data = linspace(0,10,100);
    %   obj.y_data = 1.5*obj.x_data + 2*randn(1,100);
    %
    %   % Plot raw data and the data ellipse
    %   obj.plot_data()
    %   hold on;
    %   obj.plot_ellipse();
    

    properties
        x_data  % X data
        y_data  % Y Data
        confidence = 0.95  % Confidence interval for the data ellipse
    end
    properties (SetAccess = private)
        chisquare     % chisquare distribution value for the confidence interval
        covariance    % Covariance matrix of [X, Y]
        eigenvalues   % Eigenvalues of the covariance matrix, sorted from smallest to larges
        eigenvectors  % Eigenvectors of the covariance matrix, corresponding to the eigenvalues

        ellipse_angle   % Orientation of the ellipse (in degrees, relative to horizontal)
        ellipse_center  % Coordinates of the center of the ellipse
        ellipse_rotation_matrix  % A transformation to rotate the ellipse
        minor_radius    % Minor radius of the ellipse
        major_radius    % Major radius of the ellipse
    end
    properties (Dependent = true)
        
        % The coordinates of the ellipse contour. 
        % Data is returned as a structure with "x" and "y" coordinates.
        %
        % Example:
        %     coords = obj.ellipse_coordinates()
        %     plot(coords.x, coords.y, '-')
        %
        ellipse_coordinates  
    end

    methods
        function obj = data_ellipse(x_data, y_data, confidence)
            % data_ellipse Constructs an instance of this class
            %   
            % X and Y data can be passed as optional arguments, as well as
            % the confidence interval.
            %
            % Optional arguments:
            %    x_data: x data
            %    y_data: y data
            %    confidence: confidence interval (default = 0.95).
            %
            % Example:
            %   % Initialize an object with x and Y data.
            %   obj95 = data_ellipse(x, y);
            %   % Initialize the object specifying a 99% confidence
            %   % interval.
            %   obj99 = data_ellipse(x, y; 99);
            %   % Initialize an object and assign data later.
            %   obj = data_ellipse();
            %   obj.x_data = x;
            %   obj.y_data = y;
            %

            if exist('x_data', 'var') && ~isempty(x_data)
                obj.x_data = x_data;
            end
            if exist('y_data', 'var') && ~isempty(y_data)
                obj.y_data = y_data;
            end
            if exist('confidence', 'var') && ~isempty(confidence)
                obj.confidence = confidence;
            end
        end
        
        function demo(~)
            % data_ellipse.demo() Runs a demonstration of how to use the
            % class. It uses synthetic data.
            %
            % Example
            %   obj = data_ellipse;
            %   obj.demo()
            %

            % Generate synthetic data
            N = 50;
            x = randn(N,1);
            y = 1.5 * x + randn(N,1);

            % Create object and set data
            obj = data_ellipse(x, y);
            
            % Plot the data            
            figure
            h = plot(x, y, 'k.');
            hold on; axis equal;

            % Plot the ellipse for 95%
            h(2) = obj.plot_ellipse();

            % Plot the ellipse for 99%
            obj.confidence = 0.99;
            h(3) = obj.plot_ellipse_contour();
            
            % Add linear regression for comparison
            h(4) = obj.plot_linear_regression();
            

            legend(h, 'Data', '95% data ellipse', '99% data ellipse', ...
                'Linear regression', ...
                'Location', 'NorthWest')

        end

        function h = plot_data(obj, varargin)
            % h = data_ellipse.plot_data() plots the original data
            % 
            % Optional input:
            %   Line specification. For example: '.r'.
            %
            % Returns:
            % 
            %    h: the graphics handle of the ellipse

            if isempty(varargin)
                LineSpec = {'.'};
            else
                LineSpec = varargin;
            end
            
            h = plot(obj.x_data, obj.y_data, LineSpec{:});

        end

        function [h, mdl] = plot_linear_regression(obj, varargin)
            % [h, mdl] = data_ellipse.plot_data() plots a linear regression 
            % of the data
            % 
            % Optional input:
            %   Line specification. For example: '.r'.
            %
            % Returns:
            % 
            %    h: the graphics handle of the ellipse
            %    mdl: the linear model
            %

            if isempty(varargin)
                LineSpec = {'.'};
            else
                LineSpec = varargin;
            end
            
            mdl = fit(obj.x_data, obj.y_data,'poly1');
            x_regression = [min(obj.x_data), max(obj.x_data)];
            y_regression = mdl(x_regression);
            h = plot(x_regression, y_regression, 'k--');

        end


        function h = plot_ellipse_contour(obj, varargin)
            % data_ellipse.plot_ellipse_contour() plots the contour of the
            % ellipse.
            % 
            % Optional input:
            %   Line specification. For example: '-.r'.
            %
            % Returns:
            % 
            %    h: the graphics handle of the ellipse

            if isempty(varargin)
                LineSpec = {'-'};
            else
                LineSpec = varargin;
            end
            coords = obj.ellipse_coordinates;
            h = plot(coords.x, coords.y, LineSpec{:});
        
        end

        function h = plot_ellipse(obj, varargin)
            % data_ellipse.plot_ellipse() plots the data ellipse as a
            % semi-transparent surface
            % Graphic characteristics of the ellipse can be passed as
            % optional arguments.
            % 
            % Returns;
            %     h: the graphic handle of the ellipse
            %
            % Example:
            %     obj = ellipse_data(x,y);
            %     h = obj.plot_ellipse('g', 'FaceAlpha', 0.5); 
            % 
            
            if isempty(varargin)
                LineSpec = {'g', 'FaceAlpha', 0.3, 'Edgecolor', 'none'};
            else
                LineSpec = varargin;
            end

            coords = obj.ellipse_coordinates;
            h = patch(coords.x, coords.y, LineSpec{:});

        end

        function h = plot_diameters(obj, varargin)
            % plot_diameters plots the major and minor diameters of the
            % ellipse.
            %
            % Returns:
            %    h: the graphic handles of the two diameters
            %
            
            if isempty(varargin)
                LineSpec = {'k-'};
            else
                LineSpec = varargin;
            end
            
            X0 = obj.ellipse_center(1);
            Y0 = obj.ellipse_center(2);
            
            % Get vectors in the correct orientation of the ellipse
            r = [obj.major_radius , 0 ; 0 , obj.minor_radius] * obj.ellipse_rotation_matrix;
            radius1 = [X0 , Y0 ; r(1,:)];
            radius2 = [X0 , Y0 ; r(2,:)];
            
            % Plot the diameters
            h    = plot([X0 - radius1(2,1) , X0 + radius1(2,1)] , [Y0- radius1(2,2) , Y0 + radius1(2,2)] , LineSpec{:});
            h(2) = plot([X0 - radius2(2,1) , X0 + radius2(2,1)] , [Y0- radius2(2,2) , Y0 + radius2(2,2)] , LineSpec{:});
            
        end

        function set.y_data(obj,y_data)
            % set.y_data y_data setter
            obj.y_data = y_data;
            obj.calculate();
        
        end
        
        function set.x_data(obj,x_data)
            % set.x_data x_data setter
            obj.x_data = x_data;
            obj.calculate();
        end

        function set.confidence(obj, confidence)
            if confidence >= 1 || confidence < 0
                error('Confidence should be between zero and one')
            end

            obj.confidence = confidence;
            obj.calculate();
        
        end

        function coords = get.ellipse_coordinates(obj)
        % data_ellipse.ellipse_coordinates returns the coordinates of the
        % data ellipse as a structure.
            
            % Define the ellipse in parametric form over 50 points.
            theta_grid = linspace(0,2*pi, 50);
            coords_tmp  = [obj.major_radius*cos(theta_grid); ...
                           obj.minor_radius*sin(theta_grid)]';

            % Rotate the ellipse
            coords_tmp = coords_tmp * obj.ellipse_rotation_matrix;
            coords.x = coords_tmp(:,1) + obj.ellipse_center(1);
            coords.y = coords_tmp(:,2) + obj.ellipse_center(2);

        end

    end
    methods (Access = private)
        function checkInput(obj)
            % checkInput() Checks for input consistency
            
            if length(obj.x_data) ~= length(obj.y_data)
                error('x_data and y_data should have the same length.')
            end
            
            % Check x data size and vector orientation
            sz = size(obj.x_data);
            if min(sz) ~= 1
                error('X and Y data should have a size of [Mx1]')
            end
            if sz(2) > sz(1)
                obj.x_data = obj.x_data';
            end

            % Check y data size and vector orientation
            sz = size(obj.y_data);
            if min(sz) ~= 1
                error('X and Y data should have a size of [Mx1]')
            end
            if sz(2) > sz(1)
                obj.y_data = obj.y_data';
            end

        end

        function calculate(obj)
            % data_ellipse.calculate Checks data consistency and performs
            % initial calculations
            
            % Only continue if all data is available
            if isempty(obj.x_data) || isempty(obj.y_data)
                return
            end

            % Check input
            obj.checkInput();

            % Calculate chisquare value
            obj.chisquare = sqrt(chi2inv(obj.confidence, 2));

            % Covariance matrix, eigenvalues and eigenvectors
            obj.covariance = cov([obj.x_data, obj.y_data]);
            [obj.eigenvectors, obj.eigenvalues] = eig(obj.covariance);

            % Sort eigenvalues from smallest to largest (if necessary)
            if ~issorted(diag(obj.eigenvalues))
                [obj.eigenvalues,idx] = sort(diag(obj.eigenvalues));
                obj.eigenvectors = obj.eigenvectors(:, idx);
            end

            % Calculate radii
            obj.minor_radius = obj.chisquare * sqrt(obj.eigenvalues(1,1));
            obj.major_radius = obj.chisquare * sqrt(obj.eigenvalues(2,2));

            % Calculate the angle between the x-axis and the largest eigenvector
            angle = atan2(obj.eigenvectors(2,2), obj.eigenvectors(1,2));

            % We want the ellipse to point right
            if(angle < 0)
                angle = angle + 2*pi;
            end

            % Store data
            obj.ellipse_angle = angle / pi * 180; % Transform in degrees
            obj.ellipse_rotation_matrix = [cos(angle), sin(angle); -sin(angle), cos(angle)];

            % Calculate the center of the ellipse
            obj.ellipse_center = mean([obj.x_data, obj.y_data]);

        end

    end
    
   methods(Hidden)
      % This section is only necessary to hide the inherited methods from 
      % the handle class in the documentation
         
      function lh = addlistener(varargin)
         lh = addlistener@handle(varargin{:});
      end
      function lh = listener(varargin)
         lh = listener@handle(varargin{:});
      end
      function notify(varargin)
         notify@handle(varargin{:});
      end
      function delete(varargin)
         delete@handle(varargin{:});
      end
      function Hmatch = findobj(varargin)
         Hmatch = findobj@handle(varargin{:});
      end
      function p = findprop(varargin)
         p = findprop@handle(varargin{:});
      end
      function TF = eq(varargin)
         TF = eq@handle(varargin{:});
      end
      function TF = ne(varargin)
         TF = ne@handle(varargin{:});
      end
      function TF = lt(varargin)
         TF = lt@handle(varargin{:});
      end
      function TF = le(varargin)
         TF = le@handle(varargin{:});
      end
      function TF = gt(varargin)
         TF = gt@handle(varargin{:});
      end
      function TF = ge(varargin)
         TF = ge@handle(varargin{:});
      end

   end    
end