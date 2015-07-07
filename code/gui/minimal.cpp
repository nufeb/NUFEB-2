/////////////////////////////////////////////////////////////////////////////
// Name:        minimal.cpp
// Purpose:     Minimal wxWidgets sample
// Author:      Julian Smart
// Modified by:
// Created:     04/01/98
// RCS-ID:      $Id$
// Copyright:   (c) Julian Smart
// Licence:     wxWindows licence
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <wx/grid.h>

// ----------------------------------------------------------------------------
// resources
// ----------------------------------------------------------------------------

// the application icon (under Windows and OS/2 it is in resources and even
// though we could still include the XPM here it would be unused)
#ifndef wxHAS_IMAGES_IN_RESOURCES
    #include "sample.xpm"
#endif

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

// Define a new application type, each program should derive a class from wxApp
class MyApp : public wxApp
{
public:
    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit();
};

// Define a new frame type: this is going to be our main frame
class MyFrame : public wxFrame
{
public:
    // ctor(s)
    MyFrame(const wxString& title);

    // event handlers (these functions should _not_ be virtual)
    void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    void Browse(wxCommandEvent& event);

private:
    // any class wishing to process wxWidgets events must use this macro
    wxDECLARE_EVENT_TABLE();

    wxTextCtrl *inputFile;
};

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// IDs for the controls and the menu commands
enum
{
    // menu items
    Minimal_Quit = wxID_EXIT,

    // it is important for the id corresponding to the "About" command to have
    // this standard value as otherwise it won't be handled properly under Mac
    // (where it is special and put into the "Apple" menu)
    Minimal_About = wxID_ABOUT,
    BUTTON_Browse = wxID_HIGHEST + 1
};

// ----------------------------------------------------------------------------
// event tables and other macros for wxWidgets
// ----------------------------------------------------------------------------

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them. It can be also done at run-time, but for the
// simple menu events like this the static method is much simpler.
wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(Minimal_Quit,  MyFrame::OnQuit)
    EVT_MENU(Minimal_About, MyFrame::OnAbout)
    EVT_BUTTON ( BUTTON_Browse, MyFrame::Browse )
wxEND_EVENT_TABLE()

// Create a new application object: this macro will allow wxWidgets to create
// the application object during program execution (it's better than using a
// static object for many reasons) and also implements the accessor function
// wxGetApp() which will return the reference of the right type (i.e. MyApp and
// not wxApp)
IMPLEMENT_APP(MyApp)

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// the application class
// ----------------------------------------------------------------------------

// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit()
{
    // call the base class initialization method, currently it only parses a
    // few common command-line options but it could be do more in the future
    if ( !wxApp::OnInit() )
        return false;

    // create the main application window
    MyFrame *frame = new MyFrame("NUFEB Model");

    // and show it (the frames, unlike simple controls, are not shown when
    // created initially)
    frame->Show(true);

    // success: wxApp::OnRun() will be called which will enter the main message
    // loop and the application will run. If we returned false here, the
    // application would exit immediately.
    return true;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

// frame constructor
MyFrame::MyFrame(const wxString& title)
       : wxFrame(NULL, wxID_ANY, title, wxDefaultPosition, wxSize(2000,1000))
{
    // set the frame icon
    SetIcon(wxICON(sample));

#if wxUSE_MENUS
    // create a menu bar
    wxMenu *fileMenu = new wxMenu;

    // the "About" item should be in the help menu
    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(Minimal_About, "&About\tF1", "Show about dialog");

    fileMenu->Append(Minimal_Quit, "E&xit\tAlt-X", "Quit this program");

    // now append the freshly created menu to the menu bar...
    wxMenuBar *menuBar = new wxMenuBar();
    menuBar->Append(fileMenu, "&File");
    menuBar->Append(helpMenu, "&Help");

    // ... and attach this menu bar to the frame
    SetMenuBar(menuBar);
#endif // wxUSE_MENUS

#if wxUSE_STATUSBAR
    // create a status bar just for fun (by default with 1 pane only)
    CreateStatusBar(2);
    SetStatusText("Enter parameters");
#endif // wxUSE_STATUSBAR


    new wxStaticText(
        this, -1, "All units should be specified in SI",
        wxPoint(0,0), wxDefaultSize,
        wxST_NO_AUTORESIZE);


    new wxStaticText(
        this, -1, "Boundary",
        wxPoint(0,40), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "x",
        wxPoint(85,40), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "y",
        wxPoint(145,40), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "z",
        wxPoint(205,40), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "Periodic",
        wxPoint(0,60), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "Fixed",
        wxPoint(0,80), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxRadioButton *periodicX = new wxRadioButton(this, -1, 
        "", wxPoint(80, 60), wxDefaultSize, wxRB_GROUP);
    wxRadioButton *fixedX = new wxRadioButton(this, -1, 
        "", wxPoint(80, 80));
    wxRadioButton *periodicY = new wxRadioButton(this, -1, 
        "", wxPoint(140, 60), wxDefaultSize, wxRB_GROUP);
    wxRadioButton *fixedY = new wxRadioButton(this, -1, 
        "", wxPoint(140, 80));
    wxRadioButton *periodicZ = new wxRadioButton(this, -1, 
        "", wxPoint(200, 60), wxDefaultSize, wxRB_GROUP);
    wxRadioButton *fixedZ = new wxRadioButton(this, -1, 
        "", wxPoint(200, 80));


    new wxStaticText(
        this, -1, "Input File:",
        wxPoint(0,120), wxDefaultSize,
        wxST_NO_AUTORESIZE);


    inputFile = new wxTextCtrl(this, -1, "", wxPoint(100,120), wxSize(400,-1),  
    wxTE_LEFT | wxTE_READONLY , wxDefaultValidator, wxTextCtrlNameStr);

    wxButton *inputFileBrowse = new wxButton(this, BUTTON_Browse, _T("Browse"), wxPoint(500,120), wxDefaultSize, 0);



    new wxStaticText(
        this, -1, "Neighbor Skin Distance (m):",
        wxPoint(0,160), wxDefaultSize,
        wxST_NO_AUTORESIZE);


    wxTextCtrl *neighbor = new wxTextCtrl(this, -1, "1.0e-5", wxPoint(200,160), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);



    new wxStaticText(
        this, -1, "Time Step (s):",
        wxPoint(0,200), wxDefaultSize,
        wxST_NO_AUTORESIZE);


    wxTextCtrl *timestep = new wxTextCtrl(this, -1, "1.0", wxPoint(100,200), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);


    new wxStaticText(
        this, -1, "HET",
        wxPoint(200,260), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "AOB",
        wxPoint(300,260), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "NOB",
        wxPoint(400,260), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "EPS",
        wxPoint(500,260), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "Maximum Growth Rate (1/s):",
        wxPoint(0,290), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "Decay Rate (1/s):",
        wxPoint(0,320), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "Yield Coefficient (-):",
        wxPoint(0,350), wxDefaultSize,
        wxST_NO_AUTORESIZE);


    wxTextCtrl *hetGrowthRate = new wxTextCtrl(this, -1, "6.944e-5", wxPoint(200,290), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *aobGrowthRate = new wxTextCtrl(this, -1, "3.4722e-5", wxPoint(300,290), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *nobGrowthRate = new wxTextCtrl(this, -1, "3.4722e-5", wxPoint(400,290), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *hetDecayRate = new wxTextCtrl(this, -1, "4.6296e-6", wxPoint(200,320), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *aobDecayRate = new wxTextCtrl(this, -1, "1.2731e-6", wxPoint(300,320), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *nobDecayRate = new wxTextCtrl(this, -1, "1.2731e-6", wxPoint(400,320), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *epsDecayRate = new wxTextCtrl(this, -1, "1.9676e-6", wxPoint(500,320), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *hetYieldRate = new wxTextCtrl(this, -1, "0.61", wxPoint(200,350), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *aobYieldRate = new wxTextCtrl(this, -1, "0.33", wxPoint(300,350), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *nobYieldRate = new wxTextCtrl(this, -1, "0.083", wxPoint(400,350), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);
    wxTextCtrl *epsYieldRate = new wxTextCtrl(this, -1, "0.18", wxPoint(500,350), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);



    new wxStaticText(
        this, -1, "Dissolved Oxygen Affinity (kg/m^3):",
        wxPoint(0,400), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "HET:",
        wxPoint(250,400), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *hetOxygen = new wxTextCtrl(this, -1, "8.1e-4", wxPoint(290,400), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);

    new wxStaticText(
        this, -1, "AOB:",
        wxPoint(410,400), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *aobOxygen = new wxTextCtrl(this, -1, "5e-4", wxPoint(450,400), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);

    new wxStaticText(
        this, -1, "NOB:",
        wxPoint(570,400), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *nobOxygen = new wxTextCtrl(this, -1, "6.8e-4", wxPoint(610,400), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);


    new wxStaticText(
        this, -1, "Nitrite Affinity (kg/m^3):",
        wxPoint(0,440), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "HET:",
        wxPoint(250,440), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *hetNitrite = new wxTextCtrl(this, -1, "3e-4", wxPoint(290,440), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);

    new wxStaticText(
        this, -1, "NOB:",
        wxPoint(410,440), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *nobNitrite = new wxTextCtrl(this, -1, "1.3e-3", wxPoint(450,440), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);

    new wxStaticText(
        this, -1, "HET Carbon Source Affinity (kg/m^3):",
        wxPoint(0,480), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *hetCarbon = new wxTextCtrl(this, -1, "1e-2", wxPoint(250,480), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);


    new wxStaticText(
        this, -1, "HET Nitrate Affinity (kg/m^3):",
        wxPoint(0,520), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *hetNitrate = new wxTextCtrl(this, -1, "3e-4", wxPoint(250,520), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);

    new wxStaticText(
        this, -1, "AOB Ammonia Affinity (kg/m^3):",
        wxPoint(0,560), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *aobAmmonia = new wxTextCtrl(this, -1, "1e-3", wxPoint(250,560), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);


    new wxStaticText(
        this, -1, "HET Reduction Factor (-):",
        wxPoint(0,600), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *hetReduction = new wxTextCtrl(this, -1, "0.6", wxPoint(250,600), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);


    new wxStaticText(
        this, -1, "Death Factor (-):",
        wxPoint(0,640), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *deathFactor = new wxTextCtrl(this, -1, "2.0", wxPoint(250,640), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);


    new wxStaticText(
        this, -1, "EPS:",
        wxPoint(0,680), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    new wxStaticText(
        this, -1, "Density (kg/m^3):",
        wxPoint(50,680), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *epsDensity = new wxTextCtrl(this, -1, "30", wxPoint(200,680), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);

    new wxStaticText(
        this, -1, "Extraction Ratio (-):",
        wxPoint(320,680), wxDefaultSize,
        wxST_NO_AUTORESIZE);

    wxTextCtrl *epsRatio = new wxTextCtrl(this, -1, "1.25", wxPoint(470,680), wxDefaultSize,  
    wxTE_LEFT , wxDefaultValidator, wxTextCtrlNameStr);







//     wxGrid *grid = new wxGrid( this,
//                     -1,
//                     wxPoint( 0, 240 ),
//                     wxSize( 400, 300 ) );

//     grid->CreateGrid( 100, 10 );
// // We can set the sizes of individual rows and columns
// // in pixels
// grid->SetRowSize( 0, 60 );
// grid->SetColSize( 0, 120 );
// // And set grid cell contents as strings
// grid->SetCellValue( 0, 0, "wxGrid is good" );
// // We can specify that some cells are read->only
// grid->SetCellValue( 0, 3, "This is read->only" );
// grid->SetReadOnly( 0, 3 );
// // Colours can be specified for grid cell contents
// grid->SetCellValue(3, 3, "green on grey");
// grid->SetCellTextColour(3, 3, *wxGREEN);
// grid->SetCellBackgroundColour(3, 3, *wxLIGHT_GREY);
// // We can specify the some cells will store numeric
// // values rather than strings. Here we set grid column 5
// // to hold floating point values displayed with width of 6
// // and precision of 2
// grid->SetColFormatFloat(5, 6, 2);
// grid->SetCellValue(0, 6, "3.1415");





}


// event handlers

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	delete inputFile;
    // true is to force the frame to close
    Close(true);
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    wxMessageBox(wxString::Format
                 (
                    "Welcome to %s!\n"
                    "\n"
                    "This is the minimal wxWidgets sample\n"
                    "running under %s.",
                    wxVERSION_STRING,
                    wxGetOsDescription()
                 ),
                 "About wxWidgets minimal sample",
                 wxOK | wxICON_INFORMATION,
                 this);
}

void MyFrame::Browse(wxCommandEvent& WXUNUSED(event))
{
    wxFileDialog* OpenDialog = new wxFileDialog(
		this, "Choose a file to open", wxEmptyString, wxEmptyString, 
		"Input Files (*.in)|*.in|Text Files (*.txt)|*.txt",
		wxFD_OPEN, wxDefaultPosition);
 
	// Creates a "open file" dialog with 4 file types
	if (OpenDialog->ShowModal() == wxID_OK) // if the user click "Open" instead of "Cancel"
	{
		inputFile->SetValue(OpenDialog->GetPath());
	}
 
	// Clean up after ourselves
	OpenDialog->Destroy();
}
