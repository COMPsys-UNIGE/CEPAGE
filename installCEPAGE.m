function installCEPAGE
disp(' ')
disp('********************************')
disp('* Welcome to CEPAGE Toolbox! *')
disp('********************************')
disp(' ')
disp('Checking mex compiler installation:')

warning off

myCCompiler = mex.getCompilerConfigurations('C','Selected');

if isempty(myCCompiler) || isempty(myCCompiler.Name)
    error('Unable to found C compiler for mex file')
else
    disp(['C compiler found: ',myCCompiler.Name,' compiler']);
end

myCppCompiler = mex.getCompilerConfigurations('C++','Selected');

if isempty(myCppCompiler) || isempty(myCppCompiler.Name)
    error('Unable to found C++ compiler for mex file')
else
    disp(['C++ compiler found: ',myCppCompiler.Name,' compiler']);
end
disp(' ')


useBoost = false;

answ = input('Do you want to use odeint c++ integrator [Y/n]? ','s');


if strcmpi(answ,'y') || isempty(answ)
    useBoost = true;
end

if ~useBoost
    fin = fopen('function/getCEPAGEPar.m','w');
    fprintf(fin,'function CEPAGEPar = getCEPAGEPar()\n');
    fprintf(fin,['CEPAGEPar.useBoost = false;\n']);
    fprintf(fin,['CEPAGEPar.boostDir = '''';']);
    fclose(fin);
    
else
    boostDir = input('Insert boost installation directory\n','s');
    
    if (boostDir(end) == '\' || boostDir(end) == '/') && numel(boostDir) ~= 1
        boostDir = boostDir(1:end-1);
    end
    
    
    disp(' ')
    disp('Checking odeint installation...')
    
    
    fout = fopen('try.cpp','w');
    fprintf(fout,['#include <iostream>\n',...
        '\n',...
        '#include <boost/array.hpp>\n',...
        '\n',...
        '#include <boost/numeric/odeint.hpp>\n',...
        '\n',...
        'using namespace std;\n',...
        'using namespace boost::numeric::odeint;\n',...
        '\n',...
        'const double sigma = 10.0;\n',...
        'const double R = 28.0;\n',...
        'const double b = 8.0 / 3.0;\n',...
        '\n',...
        'typedef boost::array< double , 3 > state_type;\n',...
        '\n',...
        'void lorenz( const state_type &x , state_type &dxdt , double t )\n',...
        '{\n',...
        '    dxdt[0] = sigma * ( x[1] - x[0] );\n',...
        '    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];\n',...
        '    dxdt[2] = -b * x[2] + x[0] * x[1];\n',...
        '}\n',...
        '\n',...
        'int main(int argc, char **argv)\n',...
        '{\n',...
        '    state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions\n',...
        '    integrate( lorenz , x , 0.0 , 25.0 , 0.1  );\n',...
        '}\n',...
        '\n']);
    fclose(fout);
    
    isOk = true;
    
    try
        eval(['mex -silent -c try.cpp -I"',boostDir,'"']);% -L"',boostDir,'/lib"'])
    catch
        isOk = false;
    end
    
    delete('try.*')
    
    if ~isOk
        fin = fopen('function/getCEPAGEPar.m','w');
        fprintf(fin,'function CEPAGEPar = getCEPAGEPar()\n');
        fprintf(fin,['CEPAGEPar.useBoost = false;\n']);
        fprintf(fin,['CEPAGEPar.boostDir = '''';']);
        fclose(fin);
        error('Unable to find odeint library');
    else
        
        fin = fopen('function/getCEPAGEPar.m','w');
        fprintf(fin,'function CEPAGEPar = getCEPAGEPar()\n');
        fprintf(fin,['CEPAGEPar.useBoost = true;\n']);
        boostDirToPrint = strrep(boostDir,'\','\\');
        
        fprintf(fin,['CEPAGEPar.boostDir = ''',boostDirToPrint,''';']);
        fclose(fin);
        
        disp('Odeint library found');
    end
    
    
    
    
    
end



disp(' ');
disp('Compiling mex files...')



% First compile libraries
cd c_file

mkdir bin;

infoFile = dir('src');

totalString = [];

for i=1:numel(infoFile)
    tmpFile = infoFile(i).name;
    if numel(tmpFile) > 4
        if strcmp(tmpFile(end-3:end),'.cpp')
            eval(['mex -silent  -c src/',tmpFile,' -outdir bin']);
            
            if ispc
                totalString = [totalString,' bin/',tmpFile(1:end-3),'obj'];
            elseif isunix
                totalString = [totalString,' bin/',tmpFile(1:end-3),'o'];
            else
                error('Unsopported operative system');
            end
        end
    end
    
end

if isunix
    status = system(['ar rcs libCEPAGE.a ',totalString]);
elseif ispc
    status = system(['ar rcs libCEPAGE.lib ',totalString]);
else
    error('Unsopported operative system');
end

mex -silent  -c explicit_eulero.cpp -L. -lCEPAGE;
mex -silent  -c explicit_eulero_delayed.cpp -L. -lCEPAGE;
mex -silent -c explicit_euleroEvents.cpp -L. -lCEPAGE;
mex -silent -c explicit_euleroEvents_delayed.cpp -L. -lCEPAGE;

mex -silent  -c implicit_eulero.cpp -L. -lCEPAGE;
mex -silent  -c implicit_euleroEvents.cpp -L. -lCEPAGE;

disp('eulero integrators compiled')



if useBoost
    eval(['mex -c odeint.cpp -silent -L. -lCEPAGE -I"',boostDir,'"']);
    eval(['mex -c odeint_delayed.cpp -silent -L. -lCEPAGE -I"',boostDir,'"']);
    eval(['mex -c odeintEvents.cpp -silent -L. -lCEPAGE -I"',boostDir,'"']);
    eval(['mex -c odeint_delayedEvents.cpp -silent -L. -lCEPAGE -I"',boostDir,'"']);

    disp('boost integrators compiled')
end

warning on

cd ..

disp(' ')
disp('Done.')
disp(' ')

disp(' ')
disp('Adding folders to path...')
disp(' ')

addpath([pwd,'/classes'])
disp([pwd,'/classes'])
addpath([pwd,'/c_file'])
disp([pwd,'/c_file'])
addpath(genpath([pwd,'/function']))
disp([pwd,'/function'])

disp(' ')
disp('Done.')
disp(' ')


status = savepath();
if status
    warning on
    warning('Path cannot be saved. You will need to reinstall the toolbox after restarting MATLAB');
else
    disp('Path succesfully saved')
end
disp(' ')

disp(' ')
disp('Installation complete.');

end


