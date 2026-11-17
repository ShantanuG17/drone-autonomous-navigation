function Pathfinding(mode, agentSteps)
% =========================================================================
% Pathfinding Simulation (single file, no toolboxes, no GUI)
%
% ENTRY POINT (this file):
%   Pathfinding                 original map + Cases 1-3, animated
%   Pathfinding('fast')         same, without animation pauses
%   Pathfinding('run',[a b c])  choose agent step for Cases 1-3
%   Pathfinding('retain')       Case 1 with agent already past the
%                                           obstacle -> route kept, no replan
%   Pathfinding('unreachable')  SYNTHETIC test map (not a prescribed
%                                           scenario) -> "no route" handling
%   Pathfinding('test')         automated verification (Section 5)
%   Pathfinding('gui')          interactive GUI: random or Fig. 1 map,
%                                           set/randomise S and G, find path, add
%                                           obstacles while the agent is moving
%
% GRID CONVENTION (Section 1)
%   map(r,c): 0 = free, 1 = blocked.  r = row, c = column (MATLAB indices).
%   Display:  x = c (left -> right), y = r (top -> bottom).
%   Row 1 is drawn at the TOP; this is set explicitly with YDir = 'reverse'
%   (no implicit change of origin). All cells are stored as [r c]:
%   S = [15 2]  (x,y) = (2,15);   G = [15 15]  (x,y) = (15,15).
%
% MOVEMENT MODEL (Section 2)
%   4-connected (up, down, left, right), cost 1 per move. No diagonals, so
%   corner-cutting cannot occur. Path cost = number of moves. The same rule
%   is used by the planner (bfs4) and the checker (validateRoute).
%
% ALGORITHM
%   Breadth-first search. With equal step costs BFS expands cells in order of
%   distance, so the first time G is dequeued its path has the fewest moves
%   (optimal). Dijkstra would give the same result with extra cost; A* would
%   save little on a corridor map of 256 cells.
%
% CODE SECTIONS: map representation | obstacle updates | pathfinding |
%   route validation | simulated movement | visualisation | tests
% =========================================================================
if nargin < 1, mode = 'run'; end
if nargin < 2, agentSteps = [4 4 4]; end   % moves made before the obstacle appears
cfg.plot = true; cfg.animate = true; cfg.dt = 0.08;

switch lower(mode)
    case {'run','fast'}
        cfg.animate = ~strcmpi(mode,'fast');
        runScenario('Original map', [], NaN, cfg);
        for k = 1:3
            runScenario(sprintf('Case %d',k), scenarioObstacles(k), agentSteps(k), cfg);
        end
    case 'retain'
        % Same prescribed Case 1 obstacle; agent already beyond it on row 9.
        runScenario('Case 1 (agent past obstacle)', scenarioObstacles(1), 15, cfg);
    case 'unreachable'
        % SYNTHETIC: Case 2 + Case 3 obstacles together close all three
        % corridors. Used only to test the "no route" behaviour.
        obs = unique([scenarioObstacles(2); scenarioObstacles(3)],'rows');
        runScenario('SYNTHETIC unreachable test', obs, 4, cfg);
    case 'test'
        runTests();
    case 'gui'
        runGUI();
    otherwise
        error('Unknown mode. See the header comment for valid modes.');
end
end

%% ======================= MAP REPRESENTATION =============================
function [map,S,G] = originalMap()
% Fig. 1 of the brief. Start fully blocked, then carve the free corridors.
map = ones(16);
map(2,2:15)  = 0;    % top corridor        (row 2)
map(5,2:15)  = 0;    % middle corridor     (row 5)
map(9,2:15)  = 0;    % lower corridor      (row 9)
map(2:15,2)  = 0;    % left corridor       (column 2, contains S)
map(2:15,15) = 0;    % right corridor      (column 15, contains G)
S = [15 2];          % [r c]
G = [15 15];         % [r c]
end

%% ========================= OBSTACLE UPDATES =============================
function cells = scenarioObstacles(k)
% Fig. 2: cells ([r c]) newly blocked in each independent test.
switch k
    case 1, cells = [9 8; 9 9];                 % lower corridor
    case 2, cells = [5 8; 5 9; 9 8; 9 9];       % middle + lower corridors
    case 3, cells = [2 8; 2 9; 9 8; 9 9];       % top + lower corridors
    otherwise, error('Scenario must be 1, 2 or 3');
end
end

function map = applyObstacles(map, cells)
% Returns a modified COPY. The caller always starts from originalMap(), so
% obstacles never accumulate between tests.
for i = 1:size(cells,1)
    map(cells(i,1),cells(i,2)) = 1;
end
end

%% =========================== PATHFINDING ================================
function [found,route,cost,nExpanded] = bfs4(map, startRC, goalRC)
% Unweighted BFS on a 4-connected grid. route = [r c] rows, start..goal.
% parent(v) stores the linear index of v's predecessor, which is all that
% is needed to reconstruct the route.
[nr,nc] = size(map);
found = false; route = zeros(0,2); cost = Inf; nExpanded = 0;
if map(startRC(1),startRC(2))==1 || map(goalRC(1),goalRC(2))==1, return; end
moves  = [-1 0; 1 0; 0 -1; 0 1];          % up, down, left, right
dist   = inf(nr,nc);                      % moves from start (Inf = undiscovered)
parent = zeros(nr,nc);
s = sub2ind([nr nc],startRC(1),startRC(2));
g = sub2ind([nr nc],goalRC(1),goalRC(2));
queue = zeros(nr*nc,1); head = 1; tail = 1;   % FIFO, each cell enqueued once
queue(tail) = s; tail = tail+1; dist(s) = 0;
while head < tail
    u = queue(head); head = head+1; nExpanded = nExpanded+1;
    if u == g, found = true; break; end
    [ur,uc] = ind2sub([nr nc],u);
    for k = 1:4
        vr = ur+moves(k,1); vc = uc+moves(k,2);
        if vr<1 || vr>nr || vc<1 || vc>nc, continue; end   % outside map
        if map(vr,vc)==1, continue; end                     % blocked
        v = sub2ind([nr nc],vr,vc);
        if ~isinf(dist(v)), continue; end                   % already discovered
        dist(v) = dist(u)+1; parent(v) = u;
        queue(tail) = v; tail = tail+1;
    end
end
if ~found, return; end
cost = dist(g);
idx = g;
while idx(1) ~= s, idx = [parent(idx(1)); idx]; end   % follow parents back to S
[r,c] = ind2sub([nr nc],idx);
route = [r c];
end

%% ========================= ROUTE VALIDATION =============================
function [ok,msg,cost] = validateRoute(map, route, startRC, goalRC)
% Independent check of a route against the movement model.
ok = false; cost = NaN;
[nr,nc] = size(map);
if isempty(route),                    msg = 'route is empty';          return; end
if ~isequal(route(1,:),startRC),      msg = 'does not begin at start'; return; end
if ~isequal(route(end,:),goalRC),     msg = 'does not end at G';       return; end
r = route(:,1); c = route(:,2);
if any(r<1 | r>nr | c<1 | c>nc),     msg = 'leaves the map';          return; end
b = find(map(sub2ind([nr nc],r,c))==1, 1);
if ~isempty(b), msg = sprintf('cell (r,c)=(%d,%d) is blocked',r(b),c(b)); return; end
stepLen = sum(abs(diff(route,1,1)),2);           % Manhattan length of each move
if any(stepLen ~= 1), msg = 'contains a move that is not 4-connected'; return; end
cost = sum(stepLen);                              % 1 per move
ok = true; msg = 'valid';
end

%% ===================== SCENARIO / SIMULATED MOVEMENT ====================
function res = runScenario(label, obsCells, agentStep, cfg)
% Section 3 sequence: plan -> advance -> obstacle -> validate -> replan -> go.
res = struct('label',label,'replanned',false,'reachedGoal',false, ...
             'stopCell',[],'travelled',zeros(0,2),'totalCost',NaN);
[map0,S,G] = originalMap();
fprintf('\n===== %s =====\n', label);

% 1. Initial route on the original map
[found,route,bfsCost] = bfs4(map0,S,G);
if ~found, error('No route on the original map (map definition error).'); end
[ok,msg,cost] = validateRoute(map0,route,S,G);
assert(ok && cost==bfsCost, 'Initial route failed validation: %s', msg);
fprintf('Initial route: %d moves, cost %d (validated).\n', size(route,1)-1, cost);
if cfg.plot
    fig = figure('Name',label,'NumberTitle','off','Position',[60 120 1350 470]);
end

% Original map: no obstacle, travel the whole route.
if isempty(obsCells)
    if cfg.plot
        ax = subplot(1,1,1,'Parent',fig);
        drawMap(ax,map0,[],S,G);
        plotRoute(ax,route,'b-',2.5);
        title(ax,sprintf('%s: route S->G, cost %d (validated)',label,cost));
        animateAgent(ax,route,cfg);
    end
    res.reachedGoal = true; res.stopCell = G; res.travelled = route; res.totalCost = cost;
    fprintf('Agent reached G. Total cost %d.\n', cost);
    return;
end

% 2. Advance the agent along the initial route (route(1,:) is S)
assert(agentStep>=1 && agentStep < size(route,1)-1, ...
       'agentStep must leave the agent strictly before G');
travelled = route(1:agentStep+1,:);
cur = travelled(end,:);
fprintf('Agent advanced %d moves to (r,c)=(%d,%d), (x,y)=(%d,%d).\n', ...
        agentStep, cur(1), cur(2), cur(2), cur(1));

% 3. Introduce the prescribed obstacle (on a fresh copy of the original map)
map1 = applyObstacles(map0, obsCells);
assert(map1(cur(1),cur(2))==0, 'Agent cell is blocked by the obstacle; choose another step');
fprintf('Obstacle added at (r,c): %s\n', mat2str(obsCells));

% 4. Check the REMAINING part of the planned route on the updated map
remaining = route(agentStep+1:end,:);
[stillOK,why] = validateRoute(map1,remaining,cur,G);

% 5. Keep the route, or replan from the current cell
if stillOK
    fprintf('Remaining route still valid -> replanning NOT required.\n');
    newRoute = remaining; haveRoute = true;
else
    fprintf('Remaining route INVALID (%s) -> replanning from current cell.\n', why);
    [haveRoute,newRoute,newBfsCost] = bfs4(map1,cur,G);
    if haveRoute
        [ok,msg,newCost] = validateRoute(map1,newRoute,cur,G);
        assert(ok && newCost==newBfsCost, 'Replanned route failed validation: %s', msg);
        fprintf('Replacement route: cost %d (validated).\n', newCost);
    else
        fprintf('NO ROUTE from (r,c)=(%d,%d) to G. Agent stays where it is.\n', cur(1), cur(2));
    end
    res.replanned = true;
end

% Visualisation of steps 1-6
if cfg.plot
    ax1 = subplot(1,3,1,'Parent',fig);
    drawMap(ax1,map0,[],S,G); plotRoute(ax1,route,'b-',2.5);
    plotRoute(ax1,travelled,'g-',4); plotAgent(ax1,cur);
    title(ax1,sprintf('1-2: initial route (cost %d), agent at (x,y)=(%d,%d)', ...
          size(route,1)-1,cur(2),cur(1)));

    ax2 = subplot(1,3,2,'Parent',fig);
    drawMap(ax2,map1,obsCells,S,G); plotRoute(ax2,travelled,'g-',4);
    if stillOK
        plotRoute(ax2,remaining,'b-',2.5);
        title(ax2,'3-4: obstacle added (orange) - remaining route VALID');
    else
        plotRoute(ax2,remaining,'r--',2.5);
        title(ax2,'3-4: obstacle added (orange) - remaining route INVALID');
    end
    plotAgent(ax2,cur);

    ax3 = subplot(1,3,3,'Parent',fig);
    drawMap(ax3,map1,obsCells,S,G); plotRoute(ax3,travelled,'g-',4);
    if haveRoute
        plotRoute(ax3,newRoute,'b-',2.5);
        if stillOK, t = 'retained'; else, t = 'replacement'; end
        title(ax3,sprintf('5-6: %s route, cost %d from current cell', t, size(newRoute,1)-1));
        animateAgent(ax3,newRoute,cfg);
    else
        plotAgent(ax3,cur);
        title(ax3,'5-6: NO ROUTE to G');
        text(ax3,8.5,8.5,{'NO ROUTE AVAILABLE',sprintf('from (x,y)=(%d,%d)',cur(2),cur(1))}, ...
             'Color','r','FontWeight','bold','FontSize',13, ...
             'HorizontalAlignment','center','BackgroundColor','w');
    end
end

% 6. Continue to G (or stop)
if haveRoute
    res.travelled = [travelled; newRoute(2:end,:)];
    res.reachedGoal = true; res.stopCell = G;
    res.totalCost = size(res.travelled,1)-1;
    fprintf('Agent reached G. Total moves travelled from S: %d.\n', res.totalCost);
else
    res.travelled = travelled; res.stopCell = cur;
end
end

function animateAgent(ax, route, cfg)
% Moves the agent marker one cell per step and leaves a green trail.
h = plotAgent(ax, route(1,:));
if ~cfg.animate, set(h,'XData',route(end,2),'YData',route(end,1)); return; end
for i = 2:size(route,1)
    plot(ax, route(i-1:i,2), route(i-1:i,1), 'g-', 'LineWidth', 4);
    uistack(h,'top');
    set(h,'XData',route(i,2),'YData',route(i,1));
    drawnow; pause(cfg.dt);
end
end

%% ============================ VISUALISATION =============================
function drawMap(ax, map, newCells, S, G)
% White = free, black = blocked, orange = newly blocked. x = c, y = r.
[nr,nc] = size(map);
rgb = repmat(1-map,[1 1 3]);                       % free 1 -> white, blocked -> black
for i = 1:size(newCells,1)
    rgb(newCells(i,1),newCells(i,2),:) = reshape([1 0.55 0],1,1,3);
end
image(ax,'XData',1:nc,'YData',1:nr,'CData',rgb);
set(ax,'YDir','reverse');                          % row 1 at the top (explicit)
axis(ax,'equal'); axis(ax,[0.5 nc+0.5 0.5 nr+0.5]); hold(ax,'on');
set(ax,'XTick',1:nc,'YTick',1:nr,'FontSize',7);
xlabel(ax,'x = column c'); ylabel(ax,'y = row r');
for k = 0.5:1:nc+0.5, plot(ax,[k k],[0.5 nr+0.5],'Color',[0.75 0.75 0.75]); end
for k = 0.5:1:nr+0.5, plot(ax,[0.5 nc+0.5],[k k],'Color',[0.75 0.75 0.75]); end
text(ax,S(2),S(1),'S','Color','r','FontWeight','bold','FontSize',12,'HorizontalAlignment','center');
text(ax,G(2),G(1),'G','Color',[0 0.45 0.9],'FontWeight','bold','FontSize',12,'HorizontalAlignment','center');
end

function plotRoute(ax, route, style, w)
plot(ax, route(:,2), route(:,1), style, 'LineWidth', w);    % x = c, y = r
end

function h = plotAgent(ax, rc)
h = plot(ax, rc(2), rc(1), 'o', 'MarkerSize', 11, ...
         'MarkerFaceColor', [0.1 0.7 0.1], 'MarkerEdgeColor', 'k');
end

%% ============================== TESTS ===================================
function runTests()
% Section 5 checks. Prints PASS/FAIL; errors on the first failure.
cfg.plot = false; cfg.animate = false; cfg.dt = 0;
[map0,S,G] = originalMap();

% T1: map reproduces Fig. 1 (independent text encoding of the figure)
fig1 = ['################'; '#..............#'; '#.############.#'; '#.############.#'
        '#..............#'; '#.############.#'; '#.############.#'; '#.############.#'
        '#..............#'; '#.############.#'; '#.############.#'; '#.############.#'
        '#.############.#'; '#.############.#'; '#.############.#'; '################'];
check(isequal(map0, double(fig1=='#')) && map0(S(1),S(2))==0 && map0(G(1),G(2))==0, ...
      'T1 original map matches Fig. 1; S and G free');

% T2: each modification = original + only its own cells (no accumulation)
for k = 1:3
    m = applyObstacles(map0, scenarioObstacles(k));
    check(nnz(m - map0) == size(scenarioObstacles(k),1) && all(m(:) >= map0(:)), ...
          sprintf('T2 Case %d adds only its own obstacle cells', k));
end

% T3: initial route valid and optimal (25 = 6 up + 13 across + 6 down)
[f,route,c] = bfs4(map0,S,G); [ok,~,vc] = validateRoute(map0,route,S,G);
check(f && ok && c==25 && vc==25, 'T3 initial route valid, cost 25');

% T4: validator rejects bad routes
check(~validateRoute(map0,[S; 14 3],S,[14 3]),        'T4a rejects blocked cell');
check(~validateRoute(map0,[3 2; 2 3],[3 2],[2 3]),     'T4b rejects diagonal move');
check(~validateRoute(map0,[2 2; 1 2],[2 2],[1 2]),     'T4c rejects blocked border');
check(~validateRoute(map0,[2 2; 2 3],[2 2],G),         'T4d rejects wrong end');
check(~validateRoute(map0,[2 2; 2 4],[2 2],[2 4]),     'T4e rejects 2-cell jump');

% T5: prescribed scenarios with agent before the obstacle -> replan, reach G
expected = [33 39 33];     % 4 moves + replacement (29 via row 5, 35 via row 2)
for k = 1:3
    r = runScenario(sprintf('Case %d',k), scenarioObstacles(k), 4, cfg);
    m = applyObstacles(map0, scenarioObstacles(k));
    ok = validateRoute(m, r.travelled, S, G);
    check(r.replanned && r.reachedGoal && ok && r.totalCost==expected(k), ...
          sprintf('T5 Case %d replans from current cell, valid journey, %d moves', k, expected(k)));
end

% T6: route retained when the obstacle is behind the agent
r = runScenario('Case 1 retain', scenarioObstacles(1), 15, cfg);
check(~r.replanned && r.reachedGoal && r.totalCost==25, 'T6 still-valid route retained');

% T7: unreachable goal handled (SYNTHETIC map, all corridors closed)
obs = unique([scenarioObstacles(2); scenarioObstacles(3)],'rows');
r = runScenario('Unreachable', obs, 4, cfg);
check(r.replanned && ~r.reachedGoal && isequal(r.stopCell,[11 2]), ...
      'T7 no route reported, agent not moved through obstacles');
fprintf('\nAll tests passed.\n');
end

function check(cond, name)
if cond, fprintf('PASS  %s\n', name); else, error('FAIL  %s', name); end
end

%% ================================= GUI ==================================
function runGUI()
% Interactive version of the simulation. Uses the SAME bfs4, validateRoute,
% originalMap and scenarioObstacles as the scripted scenarios.
% State is shared between the nested callbacks below.
map = []; S = []; G = [];               % occupancy grid and [r c] endpoints
route = zeros(0,2); k = 1;              % current planned route, agent = route(k,:)
agent = []; travelled = zeros(0,2);     % agent cell and cells visited so far
added = [];                             % mask: cells blocked AFTER the last plan (orange)
oldRoute = zeros(0,2);                  % last invalidated remaining route (red dashed)
running = false; stopReq = false;

fig = figure('Name','Assignment 3 - Pathfinding GUI','NumberTitle','off', ...
             'Position',[60 60 1040 680]);
ax = axes('Parent',fig,'Units','pixels','Position',[340 60 640 580]);

% ---- controls (left column) ----
lbl(fig,'Grid size:',[15 640 90 20]);
gridEdit = uicontrol(fig,'Style','edit','String','16','Position',[110 640 50 22]);
lbl(fig,'Density:',[170 640 60 20]);
densEdit = uicontrol(fig,'Style','edit','String','0.25','Position',[235 640 50 22]);
lbl(fig,'Seed (blank = random):',[15 612 150 20]);
seedEdit = uicontrol(fig,'Style','edit','String','','Position',[170 612 60 22]);
btn(fig,'Generate Random Map',[15 575 150 30],@(s,e) onGenerate());
btn(fig,'Load Fig. 1 Map',[175 575 140 30],@(s,e) onLoadFig1());

lbl(fig,'Start (x,y):',[15 540 80 20]);
startEdit = uicontrol(fig,'Style','edit','String','','Position',[95 540 60 22]);
lbl(fig,'Goal (x,y):',[170 540 70 20]);
goalEdit = uicontrol(fig,'Style','edit','String','','Position',[245 540 60 22]);
btn(fig,'Set S/G from Boxes',[15 503 150 30],@(s,e) onSetFromBoxes());
btn(fig,'Randomise S/G',[175 503 140 30],@(s,e) onRandomise());
lbl(fig,'Mouse click on map:',[15 470 130 20]);
clickMode = uicontrol(fig,'Style','popupmenu','Position',[150 470 165 24], ...
    'String',{'Toggle obstacle','Set start S','Set goal G'});

btn(fig,'Find Path (agent at S)',[15 425 300 32],@(s,e) onFindPath());

lbl(fig,'Prescribed obstacle (Fig. 1 map):',[15 390 220 20]);
casePop = uicontrol(fig,'Style','popupmenu','Position',[15 362 150 24], ...
    'String',{'Case 1','Case 2','Case 3'});
btn(fig,'Apply Case',[175 360 140 28],@(s,e) onApplyCase());
btn(fig,'Drop Random Obstacle Ahead',[15 322 300 30],@(s,e) onDropAhead());

btn(fig,'Run Simulation',[15 278 145 34],@(s,e) onRun());
btn(fig,'Stop',[170 278 145 34],@(s,e) onStop());
lbl(fig,'Speed (s per move):',[15 248 130 20]);
speed = uicontrol(fig,'Style','slider','Min',0.01,'Max',0.5,'Value',0.12, ...
    'Position',[150 248 165 20]);
logBox = uicontrol(fig,'Style','listbox','Position',[15 20 300 220], ...
    'String',{},'FontName','FixedWidth','FontSize',8);

logMsg('Generate a random map or load the Fig. 1 map.');
redraw();

% ======================= callbacks (nested) ===========================
    function onGenerate()
        if running, logMsg('Stop the simulation first.'); return; end
        n = round(str2double(get(gridEdit,'String')));
        p = str2double(get(densEdit,'String'));
        if ~isfinite(n) || n<4 || n>60, logMsg('Grid size must be 4..60.'); return; end
        if ~isfinite(p) || p<0 || p>0.9, logMsg('Density must be 0..0.9.'); return; end
        seed = str2double(get(seedEdit,'String'));
        if isfinite(seed), rng(seed); else, rng('shuffle'); end
        map = double(rand(n) < p);                 % 1 = blocked, 0 = free
        S = []; G = []; set(startEdit,'String',''); set(goalEdit,'String','');
        clearPlan();
        logMsg(sprintf('Random %dx%d map, density %.2f, %d blocked cells.', n, n, p, nnz(map)));
        logMsg('Set S and G: boxes, mouse click, or Randomise.');
        redraw();
    end

    function onLoadFig1()
        if running, logMsg('Stop the simulation first.'); return; end
        [map,S,G] = originalMap();
        showSG(); clearPlan();
        logMsg('Fig. 1 map loaded (original, no added obstacles).');
        redraw();
    end

    function onSetFromBoxes()
        if running, logMsg('Stop the simulation first.'); return; end
        if isempty(map), logMsg('Generate or load a map first.'); return; end
        [okS,s] = parseXY(get(startEdit,'String'));
        [okG,g] = parseXY(get(goalEdit,'String'));
        if ~okS || ~okG, return; end
        if map(s(1),s(2))==1 || map(g(1),g(2))==1, logMsg('S and G must be free cells.'); return; end
        if isequal(s,g), logMsg('S and G must differ.'); return; end
        S = s; G = g; clearPlan();
        logMsg(sprintf('S = (x,y)=(%d,%d), G = (%d,%d).', S(2),S(1),G(2),G(1)));
        redraw();
    end

    function onRandomise()
        if running, logMsg('Stop the simulation first.'); return; end
        if isempty(map), logMsg('Generate or load a map first.'); return; end
        [fr,fc] = find(map==0);
        if numel(fr) < 2, logMsg('Fewer than two free cells.'); return; end
        n = size(map,1); best = [];
        for t = 1:300      % prefer a connected pair at least n cells apart
            i = randperm(numel(fr),2); a = [fr(i(1)) fc(i(1))]; b = [fr(i(2)) fc(i(2))];
            if sum(abs(a-b)) < n && t < 200, continue; end
            if bfs4(map,a,b), best = [a; b]; break; end
        end
        if isempty(best)
            i = randperm(numel(fr),2); best = [fr(i(1)) fc(i(1)); fr(i(2)) fc(i(2))];
            logMsg('No connected pair found; S and G may be unreachable.');
        end
        S = best(1,:); G = best(2,:); showSG(); clearPlan();
        logMsg(sprintf('Random S = (x,y)=(%d,%d), G = (%d,%d).', S(2),S(1),G(2),G(1)));
        redraw();
    end

    function onFindPath()
        if running, logMsg('Stop the simulation first.'); return; end
        if isempty(map) || isempty(S) || isempty(G), logMsg('Need a map, S and G.'); return; end
        clearPlan(); agent = S; travelled = S;
        [f,rt,c] = bfs4(map,S,G);
        if ~f
            logMsg('NO ROUTE from S to G on this map.');
        else
            [ok,msg,vc] = validateRoute(map,rt,S,G);
            assert(ok && vc==c, 'Route failed validation: %s', msg);
            route = rt; k = 1;
            logMsg(sprintf('Route S->G: %d moves (validated).', c));
        end
        redraw();
    end

    function onApplyCase()
        if isempty(map) || ~isequal(size(map),[16 16])
            logMsg('Cases apply to the 16x16 Fig. 1 map. Load it first.'); return;
        end
        cells = scenarioObstacles(get(casePop,'Value'));
        for i = 1:size(cells,1)
            blockCell(cells(i,:));
        end
        logMsg(sprintf('Case %d obstacle applied: %s (r,c).', get(casePop,'Value'), mat2str(cells)));
        if ~running, checkAndReplan(); end
        redraw();
    end

    function onDropAhead()
        if isempty(route) || isempty(agent) || isequal(agent,G)
            logMsg('No remaining route to block.'); return;
        end
        ahead = route(k+1:end-1,:);              % never the agent cell or G
        if isempty(ahead), logMsg('Agent is next to G; nothing to block.'); return; end
        rc = ahead(randi(size(ahead,1)),:);
        blockCell(rc);
        logMsg(sprintf('Obstacle dropped on route at (x,y)=(%d,%d).', rc(2), rc(1)));
        if ~running, checkAndReplan(); end
        redraw();
    end

    function onMapClick(~,~)
        if isempty(map), return; end
        cp = get(ax,'CurrentPoint'); c = round(cp(1,1)); r = round(cp(1,2));
        if r<1 || r>size(map,1) || c<1 || c>size(map,2), return; end
        switch get(clickMode,'Value')
            case 1      % toggle obstacle (allowed while running)
                if map(r,c)==0
                    if blockCell([r c])
                        logMsg(sprintf('Obstacle added at (x,y)=(%d,%d).', c, r));
                        if ~running, checkAndReplan(); end
                    end
                else
                    map(r,c) = 0; if ~isempty(added), added(r,c) = false; end
                    logMsg(sprintf('Cell (x,y)=(%d,%d) cleared.', c, r));
                end
            case {2,3}  % set S or G (not while running)
                if running, logMsg('Stop the simulation to move S or G.'); return; end
                if map(r,c)==1, logMsg('That cell is blocked.'); return; end
                if get(clickMode,'Value')==2, S = [r c]; else, G = [r c]; end
                showSG(); clearPlan();
                logMsg(sprintf('S = (%s), G = (%s) as (x,y).', xyStr(S), xyStr(G)));
        end
        redraw();
    end

    function onRun()
        if running, return; end
        if isempty(map) || isempty(S) || isempty(G), logMsg('Need a map, S and G.'); return; end
        if isempty(agent) || isequal(agent,G), onFindPath(); end     % (re)start from S
        if isempty(route) && ~checkAndReplan(), redraw(); return; end
        running = true; stopReq = false;
        logMsg('Simulation running. Click the map to add obstacles.');
        while ~stopReq && ishghandle(fig) && ~isequal(agent,G)
            if ~checkAndReplan(), break; end      % step 4-5: validate, replan if needed
            k = k + 1; agent = route(k,:);          % step 6: one move along the route
            travelled(end+1,:) = agent; %#ok<AGROW>
            redraw(); drawnow; pause(get(speed,'Value'));
        end
        running = false;
        if ~ishghandle(fig), return; end
        if isequal(agent,G)
            logMsg(sprintf('Reached G. Moves travelled from S: %d.', size(travelled,1)-1));
        elseif stopReq
            logMsg(sprintf('Stopped at (x,y)=(%s). Run resumes from here.', xyStr(agent)));
        end
        redraw();
    end

    function onStop()
        stopReq = true;
    end

% ======================= helpers (nested) =============================
    function ok = checkAndReplan()
        % Validate the REMAINING route from the agent; replan only if invalid.
        ok = true;
        if isempty(agent) || isequal(agent,G), return; end
        if ~isempty(route)
            [valid,why] = validateRoute(map,route(k:end,:),agent,G);
            if valid, return; end
            oldRoute = route(k:end,:);
            logMsg(sprintf('Remaining route INVALID (%s).', why));
        end
        [f,rt,c] = bfs4(map,agent,G);
        if f
            [v,msg,vc] = validateRoute(map,rt,agent,G);
            assert(v && vc==c, 'Replanned route failed validation: %s', msg);
            route = rt; k = 1;
            logMsg(sprintf('Replanned from (x,y)=(%s): %d moves.', xyStr(agent), c));
        else
            route = zeros(0,2); k = 1; ok = false;
            logMsg(sprintf('NO ROUTE from (x,y)=(%s). Agent stays.', xyStr(agent)));
        end
    end

    function ok = blockCell(rc)
        % Block a cell unless it is the agent's cell, S before departure, or G.
        ok = false;
        if (~isempty(agent) && isequal(rc,agent)) || isequal(rc,G) || ...
           (isempty(agent) && isequal(rc,S))
            logMsg(sprintf('Cannot block (x,y)=(%s): agent/S/G cell.', xyStr(rc))); return;
        end
        if map(rc(1),rc(2))==1, return; end
        map(rc(1),rc(2)) = 1;
        if isempty(added), added = false(size(map)); end
        added(rc(1),rc(2)) = true; ok = true;
    end

    function clearPlan()
        route = zeros(0,2); k = 1; agent = []; travelled = zeros(0,2);
        oldRoute = zeros(0,2); added = false(size(map));
    end

    function redraw()
        if ~ishghandle(ax), return; end
        cla(ax);
        if isempty(map), title(ax,'No map'); return; end
        [nr,nc] = size(map);
        rgb = repmat(1-map,[1 1 3]);                       % white free, black blocked
        [ar,ac] = find(added & map==1);
        for i = 1:numel(ar), rgb(ar(i),ac(i),:) = reshape([1 0.55 0],1,1,3); end
        image(ax,'XData',1:nc,'YData',1:nr,'CData',rgb,'ButtonDownFcn',@onMapClick);
        set(ax,'YDir','reverse');                          % row 1 at top (explicit)
        axis(ax,'equal'); axis(ax,[0.5 nc+0.5 0.5 nr+0.5]); hold(ax,'on');
        if nc <= 40
            for t = 0.5:1:nc+0.5, line(ax,[t t],[0.5 nr+0.5],'Color',[.75 .75 .75],'HitTest','off'); end
            for t = 0.5:1:nr+0.5, line(ax,[0.5 nc+0.5],[t t],'Color',[.75 .75 .75],'HitTest','off'); end
        end
        if size(oldRoute,1) > 1
            line(ax,oldRoute(:,2),oldRoute(:,1),'Color','r','LineStyle','--','LineWidth',2,'HitTest','off');
        end
        if size(travelled,1) > 1
            line(ax,travelled(:,2),travelled(:,1),'Color',[0 .7 0],'LineWidth',4,'HitTest','off');
        end
        if size(route,1) > 1
            ahead = route(k:end,:);
            line(ax,ahead(:,2),ahead(:,1),'Color','b','LineWidth',2.5,'HitTest','off');
        end
        if ~isempty(S), text(ax,S(2),S(1),'S','Color','r','FontWeight','bold','FontSize',12, ...
                'HorizontalAlignment','center','HitTest','off'); end
        if ~isempty(G), text(ax,G(2),G(1),'G','Color',[0 .45 .9],'FontWeight','bold','FontSize',12, ...
                'HorizontalAlignment','center','HitTest','off'); end
        if ~isempty(agent)
            line(ax,agent(2),agent(1),'Marker','o','MarkerSize',11,'LineStyle','none', ...
                 'MarkerFaceColor',[.1 .7 .1],'MarkerEdgeColor','k','HitTest','off');
        end
        if isempty(route) && ~isempty(agent) && ~isequal(agent,G)
            text(ax,(nc+1)/2,(nr+1)/2,'NO ROUTE AVAILABLE','Color','r','FontWeight','bold', ...
                 'FontSize',14,'HorizontalAlignment','center','BackgroundColor','w','HitTest','off');
        end
        if isempty(route), remStr = 'none'; else, remStr = sprintf('%d moves', size(route,1)-k); end
        title(ax,sprintf('%dx%d | remaining route: %s | travelled: %d | blue=route, green=travelled, red dashed=invalidated, orange=new obstacle', ...
              nr, nc, remStr, max(size(travelled,1)-1,0)),'FontSize',8);
        xlabel(ax,'x = column c'); ylabel(ax,'y = row r');
    end

    function showSG()
        set(startEdit,'String',xyStr(S)); set(goalEdit,'String',xyStr(G));
    end

    function [ok,rc] = parseXY(str)
        % Box text "x,y" -> stored [r c] = [y x]
        ok = false; rc = [];
        v = str2double(regexp(char(str),'\s*,\s*','split'));
        if numel(v)~=2 || any(~isfinite(v)), logMsg('Use "x,y" in the boxes.'); return; end
        v = round(v);
        if any(v<1) || v(1)>size(map,2) || v(2)>size(map,1), logMsg('Point outside the map.'); return; end
        ok = true; rc = [v(2) v(1)];
    end

    function s = xyStr(rc)
        if isempty(rc), s = ''; else, s = sprintf('%d,%d', rc(2), rc(1)); end
    end

    function logMsg(str)
        fprintf('%s\n', str);
        if ~exist('logBox','var') || ~ishghandle(logBox), return; end
        L = get(logBox,'String'); if ischar(L), L = cellstr(L); end
        L{end+1} = str; if numel(L) > 300, L = L(end-299:end); end
        set(logBox,'String',L,'Value',numel(L));
    end
end

function h = lbl(fig,str,pos)
h = uicontrol(fig,'Style','text','String',str,'Position',pos,'HorizontalAlignment','left');
end

function h = btn(fig,str,pos,cb)
h = uicontrol(fig,'Style','pushbutton','String',str,'Position',pos,'Callback',cb);
end
