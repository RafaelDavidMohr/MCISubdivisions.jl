function view(c::SparseVec, inds...)
    rem = findall(i -> i in inds, c.inds)
    return view(c.cfs, inds)
end
