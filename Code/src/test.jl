using Base.Libc
i = [1,2,3,4]
j = [1,2,3,4]
for item in i
    if item > 5
        @goto escape_label
    end
end
for item in j
    if item > 5
        @goto escape_label
    end

end
println("JAAAAAAAAAA")
@label escape_label
    print("NOOO")
