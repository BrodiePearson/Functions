function removeCIbound(ax)
children = get(ax, 'children');
delete(children(1));
delete(children(2));
end