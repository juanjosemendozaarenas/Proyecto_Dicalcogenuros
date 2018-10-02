

%======================================================================
%> @ingroup matscripts
%> Creates a structure that represents a product MPS network, that can then be loaded uses tntNetworksLoad().
%> Note that the created wave function may not be normalised.
%> The variation of this function with respect to the original tntMatCreateProdMps is that here the physical
%> dimension 'd' can change from site to site. If this is indeed the case, qnums is empty
%> (for now, for simplicity), so no quantum numbers are set. 
%>
%> @param cfg A cell, each entry containing a vector representing the state on each site
%>
%> @retval wf A structure representing the product MPS
%======================================================================

function wf = tntMatCreateProdMps_InhomPhysDim_NoSymm(cfg)

L = length(cfg);

for loop = 1:L
    
    d = length(cfg{loop});

    tensor.elems_type = 'values';
    tensor.elems = cfg{loop}/sqrt(cfg{loop}'*cfg{loop});
    tensor.qn_info.qn_dir = [0,0,0];
    tensor.qn_info.qn_index = {[],[],[]};

    tensor.dims = [1,d,1];
    
    wf.nodes(loop).tensor = tensor;
    wf.nodes(loop).ids = 'LRD';
    wf.nodes(loop).indices = {0,2,1};
end

wf.start = 0;
wf.start_leg = 'L';
wf.end = L-1;
wf.end_leg = 'R';

wf.connections = cell(L,L);

for loop=2:L
    wf.connections{loop-1,loop} = 'RL';
end

end