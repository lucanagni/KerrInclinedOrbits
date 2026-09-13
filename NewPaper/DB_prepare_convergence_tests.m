function DB_prepare_convergence_tests(inputDB,grid_nx,grid_ny)
    % Prepare convergence tests parfiles to sun on server
    % inputDB: parameters to run dynamics
    % grid_nx: array of ints, radial grid sizes
    % grid_ny: array of ints, angular grid sizes

    inputDB.grid_nx = grid_nx;
    inputDB.teuk_output = 1;
    inputDB.mpi_xsize = 4;
    inputDB.mpi_ysize = 4;
    for j=1:length(grid_nx)
        inputDB.grid_nx = grid_nx(j);
        for i=1:length(grid_ny)
            inputDB.grid_ny = grid_ny(i);
            DB = DB_class(inputDB);
        end
    end

return
