function sensible_data_check(X,prompt)
if nargin<2
    prompt = '';
end
if any(isnan(X(:)))
    warning(['NaN produced, ',prompt, '. NaN count: ', num2str(sum(isnan(X(:))))]);
end
if any(isinf(X(:)))
    warning(['inf produced, ',prompt])
end
if any(X(:)==0)
    warning(['0 produced, ',prompt])
end

if any(X(:)>mean(X(:))+30*std(X(:)))
    warning(['weird outlier produced? ',prompt])
end
if any(X(:)<mean(X(:))-30*std(X(:)))
    warning(['weird outlier produced? ',prompt])
end

end