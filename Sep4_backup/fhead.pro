pro fhead, filename, extension

if n_params() eq 1 then extension=0 else extension=long(extension)

f=mrdfits(filename, extension)

help, f

end
