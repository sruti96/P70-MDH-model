function check_file_existence(url::String)

    if (isfile(url) == false)
        throw(ArgumentError("opening file $(url): No such file or directory"))
    end
end
