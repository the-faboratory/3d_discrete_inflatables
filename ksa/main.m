%Summary-------------------------------------------------------------------
% Code accompanying the paper Programming 3D Curves with Discretely
% Constrained Cylindrical Inflatables, Advanced Materials, 2023.
% Authors: Robert Baines, Sree Kalyan Patiballa, Benjamin Gorissen, Katia
% Bertoldi, and Rebecca Kramer-Bottiglio.

% What does this code do?
% Estimates positions/orientations/lengths of patches on undeformed balloon
% that cause it to assume a certain curve when inflated. 

%Inputs/Outputs------------------------------------------------------------
% In: initial balloon length and radius, file path to centerline xyz
% coordinates, and user-defined threshold values.
% See 1. USER DEFINED INPUTS 
% Out: Design pararmeters, via printing to the command line
% See 5. Report design parameters

% Idea behind Kinmatic Segmentation Algorithm------------------------------
% Given an input space curve, we look for sections
% that exceed a value of curvature, and ascribe a patch to those sections.
% Then we "deflate" the final space curve by the stretch factor determined
% by its arc length divided by the intial uninflated balloon geometry, in
% order to approximate the axial centerpoint of the patch on an
% uninflated balloon. Assumes homogenous stretches along axial
% length of balloon, which is not physically true.
% To get patch angle, we assume that the angle sweep of the binormal vector
% in the TNB frame through a section with a patch is related to the angle
% of that patch.
% Again, this is assuming that inflation doesn't change the angle of
% the patch, which is not physically true.
% To get circumfrential midpoint of each patch, we look at the angle sweep
% of the binormal vector in the TNB frame throughout the entire curve up to
% the calculated axial midpoints of each patch. The degree twist up to the
% midpoints indicates relative circumfrential midpoint in the deformed
% configuration. We then subtract half of any twist angle caused by the
% patch in question to get its relative circumfrential midpoints in the
% undeformed config. The degree of twist is calculated by either
% integrating the torsion, for non-planar segments of the curve, or for
% planar segments, using cross products and dot products to extract a
% signed angle.
%--------------------------------------------------------------------------

% Last update: 
% Robert Baines 11/18/2022
% robert.baines@yale.edu 

%% 1. USER-DEFINED INPUTS--------------------------------------------------

clear;
%close all;
l_0 =  292;                         %USER: [mm] Undeformed length of balloon
k_tol = 0.01;                       %USER: [1/mm] At this magnitude or above, the curvature is enough that we assume a patch is on that segment. 0.01 was good for all cases we tested. 
distinct_patch_spacing = 20;        %USER: [mm] Step size between two points (that may be part of a collection of close by points that represent a patch) that is big enough to indicate that there is a dinstinct new patch. Higher = fewer patches because we group more runs of points together. Lower = more patches. If too low, and neighbhoring patches are determined to lie right next to each other on the arclength, then we will get an error in KSA below.
r_bal_0 = 3.18;                     %USER: [mm] radius of undeformed balloon
path_to_centerline = 'req_centerline_test_scaled.csv'; %USER: name of centerline file (.csv format) 

%% 2. Load centerline curve------------------------------------------------

centerline = csvread(path_to_centerline);

%% 3. Calculate curvature and torsion of curve-----------------------------

% calculate arclength of centerline
[arclen_deformed,~] = arclength(centerline(:,1) ,centerline(:,2) ,centerline(:,3));
[L] = cumulative_arc_len(centerline);

% get TNB frame and geometric curvature k and torsion t
[T,N,B,k,t] = frenet(centerline(:,1),centerline(:,2),centerline(:,3));


%%%%%%
% % optional processing of input curve for smoothness: 
% gradt = gradient(t); 
% % process k and t to remove any sharp jumps
% for i = 1:length(t)
%     if abs(gradt(i)) > 0.002
%        t(i+1)=t(i); 
%     end
%     if abs(k(i)) > 100 && i ~=1
%         k(i) = k(i-1);
%     end
%     if abs(t(i)) > 100 && i ~=1
%         t(i) = t(i-1);
%     end
% end

% % note, some curves may have first and last few entries in k and t that
% % experience weird boundary effects due to derivative operations.
% % In these cases, it is suggested to pad first few and last few entries.
% % I.e.:
% k(1:6) = [k(6); k(6); k(6); k(6); k(6); k(6)];
% k((end-5):end) = [k(end-6); k(end-6); k(end-6); k(end-6); k(end-6);k(end-6)];
% t(1:6) = [t(6); t(6); t(6); t(6); t(6);t(6)];
% t((end-5):end) = [t(end-6); t(end-6); t(end-6); t(end-6); t(end-6); t(end-6)];
%%%%%%

%plot geometric torsion and curvature against arc length parameter L
figure(); hold on; yyaxis left;  plot(L,t, 'linewidth', 2); ylim([-0.1 0.1]); 
ylabel('Torsion (1/mm)'); yyaxis right; plot(L,k,'linewidth', 2);
ax = gca; ax.FontSize = 20;
ax.LineWidth = 2; box(ax,'on');
xlabel('Arc Length (mm)'); ylabel('Curvature (1/mm)')
title('Torsion and curvature vs. cumulative curve length')
% plot the curvature threshold to be considered a patch as a line
plot(L, repmat(k_tol, length(L), 1), 'k');

% Plot the 3D curve
figure();
plot_tube( centerline', 3, 50,  'join_radius', 0.1, 'face_alpha', 1); 
ax = gca; ax.FontSize = 20; camlight('headlight')
hold on;
ax = gca; ax.FontSize = 20;
ax.LineWidth = 2; xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');  grid off; axis equal;

% Add Frenet-serret frames
quiver3(centerline(:,1)',centerline(:,2)',centerline(:,3)',T(:,1)',T(:,2)',T(:,3)','color','r')
h1 = quiver3(centerline(:,1)',centerline(:,2)',centerline(:,3)',N(:,1)',N(:,2)',N(:,3)','color','g');  % normal frame corresponds to 'direction' of curvature vector
set(h1,'AutoScale','on', 'AutoScaleFactor', 0.5);
h2 = quiver3(centerline(:,1)',centerline(:,2)',centerline(:,3)',B(:,1)',B(:,2)',B(:,3)','color','b');
set(h2,'AutoScale','on', 'AutoScaleFactor', 0.5);


%% 4. Kinematic Segmentation Algorithm (KSA)----------------------------

% calculating axial midpoint---assume uniform stretch throghout the balloon
axial_stretch = arclen_deformed/l_0;

% get indicies of all points that surpasses a given magnitude of curvature
for i = 1:length(k)
    if abs(k(i)) > k_tol
        patch_bool(i) = 1;  % store 1 in this array if there is a patch.
    else
        patch_bool(i) = 0;
    end
end

% get the arc length points of the segments that contain a patch (i.e. a
% 1) in patch_bool. Result is going to be a single array with groupings of
% similar values. I.e. 40 41 42 43 44 70 71 72 73 150 151 152. This
% particular sequence would be like 3 distinct patches. Note, units in mm
patch_segments_loc = L(patch_bool==1); %Equivalent of "arr" in the algorithm block in the paper.
patch_segments_torsion = t(patch_bool == 1);

% plot representation of patch general locations atop the deformed curve
for i = 1:length(k)
    if patch_bool(i) == 1
        plot3(centerline(i,1), centerline(i,2), centerline(i,3), 'ko',  'markersize', 10)
    end
end

% Go through patch_segments_loc array and pull out all distinct runs, and
% calculate the axial midpoint of those runs as the midpoint of a patch in
% the inflated configuration.
% To do this, we compare two adjacent elements in the array and see if
% there is a reasonable step size away to be a distinct patch (i.e. 5 mm.)
% If it's greater than that step size, we assume it's another distinct
% patch section.
% Also -- we find total twist (angle rad) that each of the distinct
% segments sweeps through. This will give insight into patch angle.
% Note I've ommitted dividing by 2pi from the original formula,
% since that normalizes per period but we just want the twist angle over
% the segment.
% Also -- we find the total twist (angle in rad) up to midpoint of each
% patch. This will give insight into circumfrential midpoint of the
% patches.

patch_cnt = 1; % indexes distinct patches. equivalent of "p" in algorithm block.
distinct_cnt = 1; % indexes length coordinate of the given distinct patch. equivalent of "d" in algorithm block.

% allocate empty arrays
deformed_patch_axial_midpoint = [];
twist_angle_segment = []; %  [rad] angle that the curvature vector (same direction as normal frame in Frenet triad) sweeps through along a given segment.
twist_angle_segment_perUnitL = []; 
patch_angle_segment = [];
deformed_segment_length = [];

for j = 2:length(patch_segments_loc)
    
    if patch_segments_loc(j) - patch_segments_loc(j-1) > distinct_patch_spacing
        
        %then it's a distinct new patch, and calculate axial midpoint and twist of segment of the
        %current run:
        deformed_patch_axial_midpoint(patch_cnt) = mean(stored_distinct_run); % equivalent of "y(p)" in the algorithm block
        deformed_segment_length(patch_cnt) = stored_distinct_run(end)-stored_distinct_run(1);  % equivalent of "deformSegLength" in algorithm block
        
        %Total twist of the space curve around the segment:
        twist_angle_segment(patch_cnt) =   trapz ( stored_distinct_run, stored_distinct_run_torsion ); % equivalent of "segTwistAngle" in algorithm block
        
        patch_angle_segment(patch_cnt) = sign(twist_angle_segment(patch_cnt))* -(asin( (deformed_segment_length(patch_cnt)./axial_stretch)   / ( sqrt ( ( twist_angle_segment(patch_cnt)*r_bal_0  )^2  +  (deformed_segment_length(patch_cnt)./axial_stretch)^2) )) - pi/2);   %  Equivalent of "Theta(p)" in algorithm block.
        % note: the above can be expressed without the negative multiplier 
        % on the asin and without the -pi/2 by just using acos. 
        % The identity is just acos(x) = -asin(x)+pi/2.
        
        patch_cnt = patch_cnt + 1;
        stored_distinct_run = [];
        stored_distinct_run_torsion = [];
        distinct_cnt = 1;
        stored_distinct_run(distinct_cnt) = patch_segments_loc(j);
        distinct_cnt = distinct_cnt+1;
        
    else
        
        % then it's part of the same patch, add it to the distinct run
        stored_distinct_run(distinct_cnt) = patch_segments_loc(j);
        stored_distinct_run_torsion(distinct_cnt) = real(patch_segments_torsion(j));
        
        distinct_cnt = distinct_cnt + 1;
        
        % if we're currently processing the last element in the array, we
        % must calculate the midpoint of the stored_distinct_run  here.
        if j == length(patch_segments_loc)
            deformed_patch_axial_midpoint(patch_cnt) = mean(stored_distinct_run);
            deformed_segment_length(patch_cnt) = stored_distinct_run(end)-stored_distinct_run(1);
            twist_angle_segment(patch_cnt) =  trapz ( stored_distinct_run, stored_distinct_run_torsion );
            twist_angle_segment_perUnitL(patch_cnt) =  twist_angle_segment(patch_cnt)/(stored_distinct_run(end)-stored_distinct_run(1));
            patch_angle_segment(patch_cnt) =   sign(twist_angle_segment(patch_cnt))* -(asin( (deformed_segment_length(patch_cnt)./axial_stretch)   / ( sqrt ( ( twist_angle_segment(patch_cnt)*r_bal_0      )^2  +  (deformed_segment_length(patch_cnt)./axial_stretch)^2) ))   - pi/2);
            
        end
    end
end


%takes care of endpoints that might exceed curvature but can't
%really constitute a patch
if length(stored_distinct_run) <= 1
    patch_cnt = patch_cnt - 1;
    disp('A single point exceeds ktol and made it as a standalone stored distinct run, so we ignore it')
end


% Estimate circumfrential midpoints in deformed configuration by looking at
% rotation of binormal frame from start of curve to midpoint of each patch.
% We must keep track of this cummulative twist, because centerline twist
% can occur outside of sections that have patches (i.e. binormal frame
% rotates to change direction of curvature, implying a patch is placed in a
% different circumfrential location). Thus the torsion of the curve only
% reflects material twist when it is along along the patched segments.

noTorCnter = 1; % define counter for segments with no torsion
cum_inflect_angle = 0 ; % add up addditional inflection point angles to add to circumfrential midpoints

% find angle twist of centerline up to deformed patch midpoints.
for i = 1:patch_cnt
    
    flag = 1;
    
    angle_possible_inflec_pt = [];
    bt_patch_midpoints_torsions_seg = [];
    bt_patch_midpoint_binormalVec_seg = [];
    
    length_in_consideration = L(L < deformed_patch_axial_midpoint(i)); % consider arclength up to the deformed axial midpoint of a patch to calculate twist
    all_torsions_in_consideration = t(L < deformed_patch_axial_midpoint(i));  % get corresponding quantities at this length in consideration
    all_binormal_vecs_in_consideration = B(L < deformed_patch_axial_midpoint(i),:); 
    all_tangent_vecs_in_consideration = T(L < deformed_patch_axial_midpoint(i),:); 
    centerlinecoords_in_consideration = centerline(L < deformed_patch_axial_midpoint(i), :);
    
    % debug visualize: plot all birnormal sections in cyan to see whch
    % sections are considered as we go along loop...
    % quiver3(centerlinecoords_in_consideration(:,1)',centerlinecoords_in_consideration(:,2)',centerlinecoords_in_consideration(:,3)',all_binormal_vecs_in_consideration(:,1)',all_binormal_vecs_in_consideration(:,2)',all_binormal_vecs_in_consideration(:,3)','color','c');
    
    twist_of_centerline_up_to_ptchmidpoint(i) =  -trapz ( length_in_consideration, all_torsions_in_consideration ) ; % equivalent of "X(p)" in algorithm block.
    % NOTE: TWIST equation is negative because of right hand coordinate
    % system, where right rotations are negative (about the tangent vector
    % of the curve).
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % catch statement: if the curve contains planar sections, i.e. sections
    % with 0 torsion, then the total twist angle integration formula with trapz
    % will be 0. But in reality, we know there might be an inflection
    % point where the curvature vector swings to point to another direction. Physically this
    % is an angle change in any subsequent pathces' circumfrential midpoints.
    % Think about a sine wave: its peaks and troughts have curvature vectors pointing opposite directions, but the curve has
    % zero torsion.
    % Thus, this section checks if the curvuature is changing direction i.e. there's
    % an inflection point using the dot product of the first binormal vector
    % in the section by all subsequent ones in that section. Since they are always in
    % same plane (0 torsion), we can tell angle between them using dot
    % product. The sign of the angle can be determined using cross products.
    % Recall a dot b = norm(a)norm(b)cos(theta).
    % Recall a x b = norm(a)norm(b)sin(theta) * n  where n is unit vector
    % perpendicular to the plane containing a and b
    
    if i ~= 1 % only proceed if not the first patch's midpoint in consideration
        
        % extract torsions and binormal vecs between current and previous patch midpoints
        % extract the tangent vectors as well. And then get the middle of
        % that set of tangent vectors, which will be used to determine the
        % angle of rotation of the binormal vector later.
        if    find( L >  prev_patch_midpoint , 1 ) > find(L < deformed_patch_axial_midpoint(i), 1, 'last' )
            disp('KSA error: too small value of distinct_patch_spacing!!!')
        end
        
        bt_patch_midpoints_torsions_seg = all_torsions_in_consideration(  find( L >  prev_patch_midpoint , 1 ) : find(L < deformed_patch_axial_midpoint(i), 1, 'last' ) );
        bt_patch_midpoint_binormalVec_seg = all_binormal_vecs_in_consideration(  find( L >  prev_patch_midpoint , 1 ) : find(L < deformed_patch_axial_midpoint(i), 1, 'last' ) , : );
        bt_patch_midpoint_tangentVec_seg = all_tangent_vecs_in_consideration(  find( L >  prev_patch_midpoint , 1 ) : find(L < deformed_patch_axial_midpoint(i), 1, 'last' ), : ) ;
        midpoint_tangent_vector = bt_patch_midpoint_tangentVec_seg( round(length(bt_patch_midpoint_tangentVec_seg)./2), :   );
        bt_patch_midpoint_coordinates = centerlinecoords_in_consideration(  find( L >  prev_patch_midpoint , 1 ) : find(L < deformed_patch_axial_midpoint(i), 1, 'last' ), :  ) ;
        
        %debug visualize: again, plot vectors of considered birornmal for
        %this between patch sections, in a black color
        %quiver3(bt_patch_midpoint_coordinates(:,1)',bt_patch_midpoint_coordinates(:,2)',bt_patch_midpoint_coordinates(:,3)',bt_patch_midpoint_binormalVec_seg(:,1)',bt_patch_midpoint_binormalVec_seg(:,2)',bt_patch_midpoint_binormalVec_seg(:,3)','color','k');
        
        % see what's the angle between first binormal vec and any of the
        % subsedquent vecs:
        for q = 1:size(bt_patch_midpoint_binormalVec_seg,1)
            
            if sum(bt_patch_midpoint_binormalVec_seg(q,:)) == 0  % if binormal vec is all 0's (ie. flat line, no curvature, angle is 0 and it thus has not rotated from previous TNB frame)
                angle_possible_inflec_pt(q) = 0 ;
                
            else % otherwise, there is some curvature being indicated by the binormal vector
                
                if ((dot(bt_patch_midpoint_binormalVec_seg(1,:),bt_patch_midpoint_binormalVec_seg(q,:),2)))== 0  % if dot product is 0, birnomal vectors are orthogonal.
                    % since dot product is 0, it gives no sign, so we must
                    % take cross product and then by dotting that cross product
                    % with the tangent vector at the midpont of the
                    % line between two patches
                    % we can tell direction of rotation: https://math.stackexchange.com/questions/285346/why-does-cross-product-tell-us-about-clockwise-or-anti-clockwise-rotation
                    pos_or_neg_90degflip = sign(dot( cross(bt_patch_midpoint_binormalVec_seg(1,:),bt_patch_midpoint_binormalVec_seg(q,:)),    midpoint_tangent_vector    ));
                    
                    angle_possible_inflec_pt(q) =  pos_or_neg_90degflip .*  (pi/2) ;
                    
                else   % could be another angle , such as 180 degrees.
                    angle_possible_inflec_pt(q) =   sign((dot(bt_patch_midpoint_binormalVec_seg(1,:),bt_patch_midpoint_binormalVec_seg(q,:),2))) .* acos(  (dot(bt_patch_midpoint_binormalVec_seg(1,:),bt_patch_midpoint_binormalVec_seg(q,:),2)) / ( norm(bt_patch_midpoint_binormalVec_seg(1,:))*norm(bt_patch_midpoint_binormalVec_seg(q,:))  )   ) ;
                    %if angle = NAN, then its dot producting with itself
                    % thus angle is 0 and the quantatiy is undefined.
                end
                
            end
            
            % add to cumulative inflection-point caused twist of centerline if angle is non-0
            if sum( bt_patch_midpoints_torsions_seg  ) == 0 && ~isnan(angle_possible_inflec_pt(q)) && abs(angle_possible_inflec_pt(q)) > 0
                disp(strcat('Segment with no torsion but an inflection point of', num2str(rad2deg(angle_possible_inflec_pt(q))), 'deg detected. Adding angle to circumfrential midpoints and incrementing to next segment!!!'));
                cum_inflect_angle = cum_inflect_angle + angle_possible_inflec_pt(q); % add cummulative angle
                twist_of_centerline_up_to_ptchmidpoint(i) = twist_of_centerline_up_to_ptchmidpoint(i) + cum_inflect_angle;
                flag = 0; % indicate we added to twistofcenterline already this loop iteration
                break % breaks out of q for loop here
                
            end
            
        end
        
        % always add previous inflection point-caused angle to
        % subsequent entries if we didn't add it already
        if flag == 1
            twist_of_centerline_up_to_ptchmidpoint(i) = twist_of_centerline_up_to_ptchmidpoint(i) + cum_inflect_angle;
        end
        flag = 1;
        
    end
    
    prev_patch_midpoint = deformed_patch_axial_midpoint(i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%% 5. Report design parameters---------------------------------------------

segment_length = deformed_segment_length./axial_stretch % [mm] estimated length of patches

undeformed_patch_axial_midpoints = real((deformed_patch_axial_midpoint)./axial_stretch) % [mm] this is the esimated axial location of the patches in the uninflated configuration
%undeformed_patch_axial_midpoints_bounds = patch_segments_loc./axial_stretch  % will give you 'bounds' for midpoint of each patch

% Below, twist_of_centerline_up_to_ptchmidpoint should give approximate
% degree locations of circumfrential patch midpoints in deformed config,
% which will be the same as the undeformed config,
% assuming there is no twist across the patches (i.e. patches are not at
% angle in undeformed). To get circumfrential midpoint when there are
% angled patches, we must subtract 0.5*twist_angle_segment from the
% configuration, essentially untwisting up to the center of the patch
twist_of_centerline_up_to_ptchmidpoint = twist_of_centerline_up_to_ptchmidpoint  - (0.5).*twist_angle_segment;

% Make first patch as circumfrential midpoint 0
if twist_of_centerline_up_to_ptchmidpoint(1) < 0
    twist_of_centerline_up_to_ptchmidpoint = twist_of_centerline_up_to_ptchmidpoint + twist_of_centerline_up_to_ptchmidpoint(1);
else
    twist_of_centerline_up_to_ptchmidpoint = twist_of_centerline_up_to_ptchmidpoint - twist_of_centerline_up_to_ptchmidpoint(1);
end

% map rotations more than 360 deg back to between 0 and 360 deg.
for b = 1:length(twist_of_centerline_up_to_ptchmidpoint)
    while abs(twist_of_centerline_up_to_ptchmidpoint(b)) >= 2*pi
        twist_of_centerline_up_to_ptchmidpoint(b) = 2*pi - abs(twist_of_centerline_up_to_ptchmidpoint(b));
    end
end

circumfrential_midpoints_deg = rad2deg(twist_of_centerline_up_to_ptchmidpoint)
circumfrential_midpoints_mm =  (r_bal_0*2*pi .*  (twist_of_centerline_up_to_ptchmidpoint./(2*pi))) % this is just the equation: unknown arc length / Circumfrence =  known degree sweep / 360

% patch orientation on undeformed geometry (theta)
patch_angles_deg = rad2deg(patch_angle_segment)
