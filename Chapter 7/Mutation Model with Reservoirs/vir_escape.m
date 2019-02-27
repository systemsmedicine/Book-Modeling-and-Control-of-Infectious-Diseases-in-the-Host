function [sig,failed]=vir_escape(vtotal,failed, thresh)
    if (vtotal > thresh) failed=2; end;
    if failed sig = 2; else sig = 1; end
end