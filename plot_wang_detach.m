function plot_wang_detach
% Plot Wang-style IBM detachment results written by wang_ibm_detach.cpp
%
% Usage
%   1. Put this file next to the result folders, or set MANUAL_DIRS below.
%   2. Run plot_wang_detach in MATLAB.
%
% Expected folder contents
%   history.csv
%   summary.txt
%   state_000000.dat ...

MANUAL_DIRS = { ...
    'result_c1_soft', ...
    'result_c1_stiff' ...
    };          % leave empty {} to auto-find folders that contain history.csv
MAKE_VIDEO  = true;
FPS         = 8;
DPI         = 300;

close all;
dirs = locateRuns(MANUAL_DIRS);
if isempty(dirs)
    error(['No result folder found. Run the C++ solver first, then set ', ...
        'MANUAL_DIRS to those output directories.']);
end
out = fullfile(fileparts(dirs{1}), 'figures');
if ~isfolder(out), mkdir(out); end

plotHistory(dirs, out, DPI);
plotSnapshots(dirs, out, DPI);
if MAKE_VIDEO
    for k = 1:numel(dirs)
        makeVideo(dirs{k}, fullfile(out, [getLast(dirs{k}) '_process.mp4']), FPS);
    end
end
fprintf('Saved MATLAB figures in: %s\n', out);
end

function dirs = locateRuns(manual)
dirs = {};
if ~isempty(manual)
    for k = 1:numel(manual)
        p = char(manual{k});
        if isfile(fullfile(p, 'history.csv'))
            dirs{end+1} = p; %#ok<AGROW>
        elseif isfile(p) && endsWith(p, 'history.csv')
            dirs{end+1} = fileparts(p); %#ok<AGROW>
        end
    end
    return
end
bases = unique({pwd, fileparts(mfilename('fullpath'))}, 'stable');
found = {};
for b = 1:numel(bases)
    d = dir(fullfile(bases{b}, '**', 'history.csv'));
    for k = 1:numel(d)
        found{end+1} = d(k).folder; %#ok<AGROW>
    end
end
dirs = unique(found, 'stable');
end

function plotHistory(dirs, out, dpi)
fig = figure('Color', 'w', 'Position', [60 60 1280 520]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
cols = lines(max(numel(dirs), 7));
for k = 1:numel(dirs)
    T = readtable(fullfile(dirs{k}, 'history.csv'), 'VariableNamingRule', 'preserve');
    nm = getLast(dirs{k});
    plot(ax1, T.time_s, T.frac_detached, 'LineWidth', 1.8, 'Color', cols(k,:), 'DisplayName', nm);
    plot(ax2, T.time_s, T.n_inter_bonds, 'LineWidth', 1.8, 'Color', cols(k,:), 'DisplayName', nm);
end
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'Detached fraction');
title(ax1, '(a) Wall-anchor loss');
legend(ax1, 'Location', 'best');
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'Living inter-unit bonds');
title(ax2, '(b) Connection fracture');
legend(ax2, 'Location', 'best');
sgtitle('Wang IBM detachment: strain-threshold fracture');
exportgraphics(fig, fullfile(out, 'wang_detach_history.png'), 'Resolution', dpi);
close(fig);
end

function plotSnapshots(dirs, out, dpi)
n = numel(dirs);
fig = figure('Color', 'w', 'Position', [40 40 420*n 720]);
tiledlayout(2, n, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n
    files = listStates(dirs{k});
    if isempty(files), continue; end
    S0 = readState(files{1});
    S1 = readState(files{end});
    ax = nexttile(k); drawState(ax, S0); title(ax, [getLast(dirs{k}) '  initial']);
    ax = nexttile(n+k); drawState(ax, S1); title(ax, [getLast(dirs{k}) '  final']);
end
sgtitle('Units (dots) and living bonds (lines), true scale');
exportgraphics(fig, fullfile(out, 'wang_detach_snapshots.png'), 'Resolution', dpi);
close(fig);
end

function makeVideo(caseDir, outFile, fps)
files = listStates(caseDir);
if numel(files) < 2
    warning('Not enough frames in %s', caseDir);
    return
end
T = [];
if isfile(fullfile(caseDir, 'history.csv'))
    T = readtable(fullfile(caseDir, 'history.csv'), 'VariableNamingRule', 'preserve');
end
v = VideoWriter(outFile, 'MPEG-4');
v.FrameRate = fps;
open(v);
guard = onCleanup(@()close(v));
fig = figure('Color', 'w', 'Position', [80 80 900 420]);
for k = 1:numel(files)
    S = readState(files{k});
    clf(fig);
    ax = axes(fig); %#ok<LAXES>
    drawState(ax, S);
    ttl = sprintf('%s   frame %d / %d', getLast(caseDir), k, numel(files));
    if ~isempty(T) && k <= height(T)
        ttl = sprintf('%s, t = %.3g s, detached = %.2f', ...
            getLast(caseDir), T.time_s(k), T.frac_detached(k));
    end
    title(ax, ttl);
    drawnow;
    writeVideo(v, getframe(fig));
end
close(fig);
clear guard
end

function files = listStates(caseDir)
d = dir(fullfile(caseDir, 'state_*.dat'));
if isempty(d), files = {}; return; end
[~, idx] = sort({d.name});
d = d(idx);
files = fullfile({d.folder}, {d.name});
end

function S = readState(fname)
% Parse the three POINT zones written by wang_ibm_detach.cpp
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open %s', fname); end
cleaner = onCleanup(@()fclose(fid));
S.field = [];
S.units = zeros(0, 2);
S.bonds = zeros(0, 4);
zone = '';
I = 0; J = 1;
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if contains(line, 'ZONE')
        zone = lower(string(line));
        It = regexp(line, 'I\s*=\s*(\d+)', 'tokens', 'once');
        Jt = regexp(line, 'J\s*=\s*(\d+)', 'tokens', 'once');
        I = str2double(It{1});
        if isempty(Jt), J = 1; else, J = str2double(Jt{1}); end
        n = I * J;
        A = textscan(fid, '%f%f%f%f%f', n, 'CollectOutput', true);
        if isempty(A) || isempty(A{1}), continue; end
        P = A{1};
        if contains(zone, 'field')
            S.field.x = reshape(P(:,1), I, J);
            S.field.y = reshape(P(:,2), I, J);
            S.field.u = reshape(P(:,3), I, J);
            S.field.v = reshape(P(:,4), I, J);
            S.field.phase = reshape(P(:,5), I, J);
        elseif contains(zone, 'unit')
            S.units = P(:,1:2);
        elseif contains(zone, 'bond')
            xy = P(:,1:2);
            nb = floor(size(xy,1) / 2);
            B = zeros(nb, 4);
            for i = 1:nb
                B(i,:) = [xy(2*i-1,1) xy(2*i-1,2) xy(2*i,1) xy(2*i,2)];
            end
            S.bonds = B;
        end
    end
end
end

function drawState(ax, S)
hold(ax, 'on');
if isstruct(S.field) && isfield(S.field, 'u')
    umag = sqrt(S.field.u.^2 + S.field.v.^2);
    pcolor(ax, S.field.x, S.field.y, umag);
    shading(ax, 'interp');
    colormap(ax, parula);
end
if ~isempty(S.bonds)
    for i = 1:size(S.bonds,1)
        plot(ax, [S.bonds(i,1) S.bonds(i,3)], [S.bonds(i,2) S.bonds(i,4)], ...
            'k-', 'LineWidth', 0.8);
    end
end
if ~isempty(S.units)
    plot(ax, S.units(:,1), S.units(:,2), 'o', ...
        'MarkerFaceColor', [0.85 0.25 0.15], 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 5);
end
axis(ax, 'equal');
axis(ax, 'tight');
box(ax, 'on');
xlabel(ax, 'x (\mum)');
ylabel(ax, 'y (\mum)');
cb = colorbar(ax);
cb.Label.String = '|u| (m/s)';
end

function name = getLast(p)
[~, name] = fileparts(char(p));
end
