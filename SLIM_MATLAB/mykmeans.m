%Standard kmeans

function group = mykmeans(V,nhat)

	dim = size(V,2);
	mu = zeros(nhat,dim);

	N = size(V,1);

	%Random initial assignments
	group = randi(nhat,N,1);

	tol = 1;
	diff = tol+1;

	while diff > tol
		old_group = group;

		% Reassign the centers
		for i=1:nhat
			mu(i,:) = mean(V(group==i,:));
		end

		%Redo the group assignments
		for i=1:N
			temp = zeros(nhat,1);
			for j=1:nhat
				temp(j) = norm(V(i,:) - mu(j,:));
			end
			[a,b] = min(temp);
			group(i) = b;
		end

		diff = nnz(group-old_group);

	end
end