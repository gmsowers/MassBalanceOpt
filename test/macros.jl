@testset "Macro tests" begin
let
    m = Model(Ipopt.Optimizer)
    @variables(m, begin
        x
        y
        z
        c1
        c2
    end)
    @set x = 1.0
    @test start_value(x) ≈ 1.0

    @set +x
    @set x = 2.0
    @test is_fixed(x)
    @test fix_value(x) ≈ 2.0

    @set -x
    @test !is_fixed(x)

    @set -x
    @test !is_fixed(x)

    @specs begin
        +c1
        +c2 
    end
    eq = @set c1 = c2
    @test is_fixed(c1)
    @test !is_fixed(c2)
    @test eq isa ConstraintRef

    @specs begin
        -x
        +y
        x ~ y
    end
    @test is_fixed(x) && !is_fixed(y)

    @bounds begin
        -x
        -y
        -z
        x > 1.0
        y < 2.0
        1.0 < z < 2.0
    end
    @test lower_bound(x) ≈ 1.0
    @test upper_bound(y) ≈ 2.0
    @test lower_bound(z) ≈ 1.0
    @test upper_bound(z) ≈ 2.0

    @bounds begin
        -Inf < z < Inf
        x > -Inf
        y < Inf
    end
    @test !has_lower_bound(z) && !has_upper_bound(z)
    @test !has_lower_bound(x)
    @test !has_upper_bound(y)

    coef = @stoic 2A + 3B => C
    @test coef[1][:A] == -2 && coef[1][:B] == -3 && coef[1][:C] == 1

    coef = @stoic begin
        2A + 3B => C
        A + 2B => 3C + 4D
    end
    @test coef[1][:A] == -2 && coef[1][:B] == -3 && coef[1][:C] == 1
    @test coef[2][:A] == -1 && coef[2][:B] == -2 && coef[2][:C] == 3 && coef[2][:D] == 4

end
end
nothing