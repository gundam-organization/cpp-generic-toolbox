//
// Created by Adrien Blanchet on 07/05/2023.
//

#ifndef SIMPLEMODMANAGER_GENERICTOOLBOX_SWITCH_BOREALIS_IMPL_H
#define SIMPLEMODMANAGER_GENERICTOOLBOX_SWITCH_BOREALIS_IMPL_H

#ifdef __SWITCH__


namespace GenericToolbox::Switch::Borealis{

  ProgressBarMonitorView::~ProgressBarMonitorView(){
    if( _execOnDelete_ ){ _execOnDelete_(); }
    brls::View::~View();
  }

  void ProgressBarMonitorView::setHeaderTitle(const std::string &header) {
    _header_ = header;
  }
  void ProgressBarMonitorView::setProgressColor(const NVGcolor &progressColor_) {
    _progressColor_ = progressColor_;
  }
  void ProgressBarMonitorView::setExecOnDelete(const std::function<void()> &execOnDelete_) {
    _execOnDelete_ = execOnDelete_;
  }

  void ProgressBarMonitorView::setTitlePtr(const std::string *titlePtr) {
    _titlePtr_ = titlePtr;
  }
  void ProgressBarMonitorView::setSubTitlePtr(const std::string *subTitlePtr) {
    _subTitlePtr_ = subTitlePtr;
  }
  void ProgressBarMonitorView::setProgressFractionPtr(double *progressFractionPtr) {
    _progressFractionPtr_ = progressFractionPtr;
  }
  void ProgressBarMonitorView::setSubProgressFractionPtr(double *subProgressFractionPtr) {
    _subProgressFractionPtr_ = subProgressFractionPtr;
  }

  void ProgressBarMonitorView::resetMonitorAddresses(){
    _titlePtr_ = nullptr;
    _subTitlePtr_ = nullptr;
    _progressFractionPtr_ = nullptr;
    _subProgressFractionPtr_ = nullptr;
  }

  void ProgressBarMonitorView::draw(
      NVGcontext *vg, int x, int y, unsigned int width, unsigned int height,
      brls::Style *style, brls::FrameContext *ctx) {

    float y_offset = 0;
    if(not _header_.empty()){
      y_offset += 12;
      // Drawing header
      nvgBeginPath(vg);
      nvgFontFaceId(vg, ctx->fontStash->regular);
      nvgFontSize(vg, style->Header.fontSize);
      nvgFillColor(vg, a(ctx->theme->textColor));
      nvgTextAlign(vg, NVG_ALIGN_LEFT | NVG_ALIGN_MIDDLE);
      nvgText(vg, x + style->Header.rectangleWidth + style->Header.padding, y - height/2., (_header_).c_str(), nullptr);

      // Separator
      nvgBeginPath(vg);
      nvgFillColor(vg, a(ctx->theme->separatorColor));
      nvgRect(vg, x, y - height/2. + style->Header.fontSize, width, 1);
      nvgFill(vg);
    }

    // Draw Title
    nvgFillColor(vg, a(ctx->theme->textColor));
    nvgFontSize(vg, style->Label.dialogFontSize);
    nvgFontFaceId(vg, ctx->fontStash->regular);
    nvgTextLineHeight(vg, 1.0f);
    nvgTextAlign(vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
    nvgBeginPath(vg);
    if(_titlePtr_ == nullptr) nvgText(vg, x + width / 2, y + y_offset + height / 2 - 1.8*style->Label.dialogFontSize, "", nullptr);
    else{ nvgText(vg, x + width / 2, y + y_offset + height / 2 - 1.8*style->Label.dialogFontSize, (*_titlePtr_).c_str(), nullptr); }

    // Draw subTitle
    if(_subTitlePtr_ != nullptr and not _subTitlePtr_->empty() ){ // .empty() does not work ?
      nvgBeginPath(vg);
      nvgFontFaceId(vg, ctx->fontStash->regular);
      nvgFontSize(vg, style->Header.fontSize);
      nvgFillColor(vg, a(ctx->theme->textColor));
      nvgTextAlign(vg, NVG_ALIGN_LEFT | NVG_ALIGN_MIDDLE);
      nvgText(vg, x + style->Header.rectangleWidth + style->Header.padding, y + y_offset + height / 2., (*_subTitlePtr_).c_str(), nullptr);
    }

    unsigned pBarYoffset = y + y_offset + 0.8f * height;
    if(_progressFractionPtr_ != nullptr){
      ProgressBarMonitorView::drawProgressBar(vg, x, width, pBarYoffset, *_progressFractionPtr_);
    }
    if( _subProgressFractionPtr_ != nullptr ){
      ProgressBarMonitorView::drawProgressBar(vg, x, width, pBarYoffset + 2 * 16, *_subProgressFractionPtr_);
    }

  }

  void ProgressBarMonitorView::drawProgressBar(NVGcontext *vg, int x, unsigned int width, unsigned int yPosition_, double fraction_) {
    if     ( fraction_ > 1 ) fraction_ = 1;
    else if( fraction_ < 0 ) fraction_ = 0;

    unsigned x_margin = 15;
    unsigned totalBarLength = width - 2*x_margin;
    unsigned progress_x_offset = x + x_margin;

    if( fraction_ != 1 ){
      /* Progress bar background */
      /* In official software this is more transparent instead of brighter. */
      nvgFillColor(vg, a(nvgRGBAf(1.f, 1.f, 1.f, 0.5f)));
      nvgBeginPath(vg);
      nvgRoundedRect(vg, float(progress_x_offset), float(yPosition_), float(totalBarLength), 16, 8);
      nvgFill(vg);
    }

    if( fraction_ != 0 ){
      /* Progress bar */
      nvgFillColor(vg, a(_progressColor_));
      nvgBeginPath(vg);
      nvgRoundedRect(vg, float(progress_x_offset), float(yPosition_), totalBarLength * fraction_, 16, 8);
      nvgFill(vg);
    }

  }


  void PopupLoadingBox::pushView(){
    // memory will be handled by brls
    _monitorView_ = new ProgressBarMonitorView();

    // creating a box will make the loading popup resized
    _loadingBox_ = new brls::Dialog(_monitorView_ );

    // make sure the user don't cancel while unfinished
    _loadingBox_->setCancelable( true );

    while( brls::Application::hasViewDisappearing() ){
      // wait for one extra frame before push
      std::this_thread::sleep_for(std::chrono::milliseconds( 16 ));
    }

    // push the box to the view
    brls::Application::pushView( _loadingBox_ );
  }
  void PopupLoadingBox::popView(){
    // is another view is disappearing? -> let it's time to go
    while( brls::Application::hasViewDisappearing() ){
      // wait for one extra frame before push
      std::this_thread::sleep_for(std::chrono::milliseconds( 16 ));
    }

    // call if it's still on top
    if( _loadingBox_ == brls::Application::getTopStackView() ){
      brls::Application::popView( brls::ViewAnimation::FADE );
    }
  }
  bool PopupLoadingBox::isOnTopView() const{
    return ( _loadingBox_ == brls::Application::getTopStackView() );
  }

  ProgressBarMonitorView *PopupLoadingBox::getMonitorView() const {
    return _monitorView_;
  }
  brls::Dialog *PopupLoadingBox::getLoadingBox() const {
    return _loadingBox_;
  }


}

#endif

#endif //SIMPLEMODMANAGER_GENERICTOOLBOX_SWITCH_BOREALIS_IMPL_H
