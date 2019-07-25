function B = converte_csr(A, IA, JA, nl)

  B = zeros(nl,nl);

  for i = 1: nl;
    p1 = IA(i)+1;
    p2 = IA(i+1);
    for p = p1:p2;
      j = JA(p)+1;
      B(i,j) = A(p);
    end 
  end
  
end     
      