function test_whiteness_test ()
    X = rand(1,1000);
    Y = filter(1,[1 -0.5],X);
    
    h = whiteness_test(X, 10);
    assert (h == 0);
    
    h = whiteness_test(Y, 10);
    assert (h == 1);
end