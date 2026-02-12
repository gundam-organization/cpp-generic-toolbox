//
// Created by Nadrino on 23/12/2023.
//

#ifndef CPP_GENERIC_TOOLBOX_VECTOR_H
#define CPP_GENERIC_TOOLBOX_VECTOR_H

// ***************************
//! Vector related tools
// ***************************

#include "GenericToolbox.Macro.h"
#include "GenericToolbox.Log.h"

#include <functional>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <cmath>
#include <list>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"


// Declaration section
namespace GenericToolbox{

  // Content management
  template<typename Elem, typename Cont> static bool isIn( const Elem& element_, const Cont& container_ );
  template<typename Elem, typename Cont, typename Lambda> static bool isIn( const Elem& element_, const Cont& container_, const Lambda& fetchElem_ );
  template<typename Elm, typename Val, typename Lambda> static int isIn(const Val& value_, const std::vector<Elm>& vector_, const Lambda& fetchElmValueFct_);

  template<typename T> static int findElementIndex( const T& element_, const std::vector<T>& vector_ );
  static int findElementIndex(const char* element_, const std::vector<std::string>& vector_ );
  template<typename T> static void insertInVector(std::vector<T> &vector_, const std::vector<T> &vectorToInsert_, size_t insertBeforeThisIndex_);
  template<typename T> static void insertInVector(std::vector<T> &vector_, const T &elementToInsert_, size_t insertBeforeThisIndex_);
  template<typename T> static void mergeInVector(std::vector<T> &vector_, const std::vector<T> &other_, bool allowDuplicates_ = true);
  template<typename T> static void addIfNotInVector(const T& element_, std::vector<T> &vector_);
  static void addIfNotInVector(const char* element_, std::vector<std::string> &vector_);
  template <typename Elm, typename Val, typename Lambda> static int findElementIndex(const Val& value_, const std::vector<Elm>& vector_, const Lambda& fetchElmValueFct_);

  // Generators
  template<typename T> static std::vector<size_t> indices(const std::vector<T> &vector_);
  template<typename T> static std::vector<T> getSubVector( const std::vector<T>& vector_, size_t beginIndex_, int endIndex_ = -1 );
  template<typename T, typename TT> static std::vector<TT> convertVectorType( const std::vector<T>& vector_, std::function<TT(T)>& convertTypeFunction_ );

  // Stats
  template <typename T> static double getAverage(const std::vector<T>& vector_, const std::function<double(const T&)>& evalElementFct_ = [](const T& var){return var;});
  template<typename T> static double getAveragedSlope(const std::vector<T> &yValues_);
  template<typename T, typename TT> static double getAveragedSlope(const std::vector<T> &yValues_, const std::vector<TT> &xValues_);
  template <typename T> static double getStdDev(const std::vector<T>& vector_, const std::function<double(const T&)>& evalElementFct_ = [](const T& var){return var;});

  // Sorting
  template <typename T, typename Lambda> static std::vector<size_t> getSortPermutation(const std::vector<T>& vectorToSort_, const Lambda& firstArgGoesFirstFct_ );
  template <typename T> static std::vector<T> getSortedVector(const std::vector<T>& unsortedVector_, const std::vector<std::size_t>& sortPermutation_);
  template <typename T> static void applyPermutation(std::vector<T>& vectorToPermute_, const std::vector<std::size_t>& sortPermutation_);
  template <typename T, typename Lambda> static void sortVector(std::vector<T>& vectorToSort_, const Lambda& firstArgGoesFirstFct_);
  template <typename T, typename Lambda> static void removeEntryIf(std::vector<T>& vector_, const Lambda& removeIfFct_);

  // Others
  template<typename T, typename TT> static T& getListEntry(std::list<T>& list_, TT index_);
  template<typename T, typename TT> static const T& getListEntry(const std::list<T>& list_, TT index_);

  // deprecated
  template <typename T> static bool doesElementIsInVector( const T& element_, const std::vector<T>& vector_ );
  static bool doesElementIsInVector(const char* element_, const std::vector<std::string>& vector_);
  template <typename Elm, typename Val, typename Lambda> static int doesElementIsInVector(const Val& value_, const std::vector<Elm>& vector_, const Lambda& fetchElmValueFct_);

  // class
  template <typename T, typename OffsetT = uint32_t>
  class CSRVector {

    static_assert(std::is_integral<OffsetT>::value, "OffsetT must be integral");
    static_assert(std::is_unsigned<OffsetT>::value, "OffsetT must be unsigned");

  public:
    using value_type  = T;
    using offset_type = OffsetT;

    class RowView;
    class RowBuilder;

    CSRVector();

    void reserve(size_t nEntries);
    void reserve_data(size_t nTotalValues);
    void clear();

    [[nodiscard]] size_t size() const;
    [[nodiscard]] bool empty() const;

    RowBuilder emplace_back();
    void push_back(const std::vector<value_type>& row);

    RowView operator[](size_t iEntry) const;

    const std::vector<offset_type>& offsets() const;
    const std::vector<value_type>& data() const;

    void append(const CSRVector& other);
    void append_move(CSRVector&& other);
    static CSRVector concat(const std::vector<CSRVector>& parts);

  private:
    std::vector<offset_type> offsets_;
    std::vector<value_type> data_;
    bool openRow_ = false;
  };

  template <typename T, typename OffsetT>
  class CSRVector<T, OffsetT>::RowView {
  public:
    RowView(const value_type* ptr, size_t n);

    const value_type& operator[](size_t k) const;

    const value_type* data() const;
    [[nodiscard]] size_t size() const;
    [[nodiscard]] bool empty() const;

    const value_type* begin() const;
    const value_type* end() const;

  private:
    const value_type* ptr_;
    size_t n_;
  };

  template <typename T, typename OffsetT>
  class CSRVector<T, OffsetT>::RowBuilder {
  public:
    RowBuilder(CSRVector* owner, size_t rowIndex);
    RowBuilder(const RowBuilder&) = delete;
    RowBuilder& operator=(const RowBuilder&) = delete;
    RowBuilder(RowBuilder&& other) noexcept;
    RowBuilder& operator=(RowBuilder&& other) noexcept;
    ~RowBuilder();

    template <class... Args>
    value_type& emplace_back(Args&&... args);

    void close();

  private:
    void close_if_needed();

    CSRVector* owner_;
    size_t rowIndex_;
    bool open_ = true;
  };


}


// Implementation section
namespace GenericToolbox {

  // Content management
  template<typename Elem, typename Cont> static bool isIn( const Elem& element_, const Cont& container_ ){
    return std::find(container_.cbegin(), container_.cend(), element_) != container_.cend();
  }
  template<typename Elem, typename Cont, typename Lambda> static bool isIn( const Elem& element_, const Cont& container_, const Lambda& fetchElem_ ){
    return std::find_if( container_.begin(), container_.end(), [&](const Elem& t){ return fetchElem_(t) == element_; }) != container_.end();
  }
  template <typename Elm, typename Val, typename Lambda> static int isIn(const Val& value_, const std::vector<Elm>& vector_, const Lambda& fetchElmValueFct_){
    return std::find_if( vector_.begin(), vector_.end(), [&](const Elm& t){ return fetchElmValueFct_(t) == value_; }) != vector_.end();
  }

  template <typename T> static bool doesElementIsInVector(const T& element_, const std::vector<T>& vector_){
    return std::find(vector_.cbegin(), vector_.cend(), element_) != vector_.cend();
  }
  static bool doesElementIsInVector(const char* element_, const std::vector<std::string>& vector_){
    return std::find(vector_.cbegin(), vector_.cend(), element_) != vector_.cend();
  }
  template <typename Elm, typename Val, typename Lambda> static int doesElementIsInVector(const Val& value_, const std::vector<Elm>& vector_, const Lambda& fetchElmValueFct_){
    return std::find_if( vector_.begin(), vector_.end(), [&](const Elm& t){ return fetchElmValueFct_(t) == value_; }) != vector_.end();
  }
  template <typename T> inline static int findElementIndex(const T& element_, const std::vector<T>& vector_ ){ // test
    auto it = std::find(vector_.cbegin(), vector_.cend(), element_);
    if( it == vector_.cend() ){ return -1; }
    return static_cast<int>( std::distance(vector_.cbegin(), it) );
  }
  inline static int findElementIndex(const char* element_, const std::vector<std::string>& vector_ ){
    auto it = std::find(vector_.cbegin(), vector_.cend(), element_);
    if( it == vector_.cend() ){ return -1; }
    return static_cast<int>( std::distance(vector_.cbegin(), it) );
  }
  template<typename T> inline static void insertInVector(std::vector<T> &vector_, const std::vector<T> &vectorToInsert_, size_t insertBeforeThisIndex_){
    if( insertBeforeThisIndex_ > vector_.size() ){
      throw std::runtime_error("GenericToolBox::insertInVector error: insertBeforeThisIndex_ >= vector_.size()");
    }
    if( vector_.empty() ){ vector_ = vectorToInsert_; return; }
    if( vectorToInsert_.empty() ){ return; }
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ <= 4)
    vector_.insert( vector_.begin() + insertBeforeThisIndex_, vectorToInsert_.begin(), vectorToInsert_.end() );
#else
    vector_.insert( vector_.cbegin() + insertBeforeThisIndex_, vectorToInsert_.cbegin(), vectorToInsert_.cend() );
#endif
  }
  template<typename T> inline static void insertInVector(std::vector<T> &vector_, const T &elementToInsert_, size_t insertBeforeThisIndex_){
    insertInVector(vector_, std::vector<T>{elementToInsert_}, insertBeforeThisIndex_);
  }
  template<typename T> inline static void mergeInVector( std::vector<T> &vector_, const std::vector<T> &other_, bool allowDuplicates_ ){
    if( allowDuplicates_ ){ GenericToolbox::insertInVector(vector_, other_, vector_.size()); }
    else{
      vector_.reserve( vector_.size() + other_.size() );
      for( auto& element : other_ ){ GenericToolbox::addIfNotInVector(element, vector_); }
      vector_.shrink_to_fit();
    }
  }
  template<typename T> static void addIfNotInVector(const T& element_, std::vector<T> &vector_){
    if( not GenericToolbox::doesElementIsInVector(element_, vector_) ){
      vector_.emplace_back(element_);
    }
  }
  static void addIfNotInVector(const char* element_, std::vector<std::string> &vector_){
    if( not GenericToolbox::doesElementIsInVector(element_, vector_) ){
      vector_.emplace_back(element_);
    }
  }
  template <typename Elm, typename Val, typename Lambda> inline static int findElementIndex(const Val& value_, const std::vector<Elm>& vector_, const Lambda& fetchElmValueFct_){
    auto it = std::find_if( vector_.begin(), vector_.end(), [&](const Elm& t){ return fetchElmValueFct_(t) == value_; });
    if( it == vector_.end() ){ return -1; }
    else{ return std::distance( vector_.begin(), it ); }
  }

  // Generators
  template<typename T> static std::vector<size_t> indices(const std::vector<T> &vector_){
    std::vector<size_t> output(vector_.size(), 0);
    for( size_t iIndex = 0 ; iIndex < output.size() ; iIndex++ ){
      output.at(iIndex) = iIndex;
    }
    return output;
  }
  template <typename T> static std::vector<T> getSubVector( const std::vector<T>& vector_, size_t beginIndex_, int endIndex_ ){
    if( endIndex_ < 0 ){ endIndex_ += vector_.size(); }
    if( beginIndex_ >= endIndex_ ){ return std::vector<T> (); }
    return std::vector<T> ( &vector_[beginIndex_] , &vector_[endIndex_+1] );
  }
  template <typename T, typename TT> static std::vector<TT> convertVectorType( const std::vector<T>& vector_, std::function<TT(T)>& convertTypeFunction_ ){
    std::vector<TT> outVec;
    for(const auto& element : vector_){
      outVec.emplace_back(convertTypeFunction_(element));
    }
    return outVec;
  }

  // Stats
  template <typename T> static double getAverage(const std::vector<T>& vector_, const std::function<double(const T&)>& evalElementFct_){
    double outVal = 0;
    for( auto& element : vector_ ){ outVal += static_cast<double>(evalElementFct_(element)); }
    return outVal / vector_.size();
  }
  template<typename T> static double getAveragedSlope(const std::vector<T> &yValues_){
    auto xValues = yValues_;
    for( size_t iVal = 0 ; iVal < yValues_.size() ; iVal++ ){
      xValues.at(iVal) = iVal;
    }
    return getAveragedSlope(yValues_, xValues);
  }
  template<typename T, typename TT> static double getAveragedSlope(const std::vector<T> &yValues_, const std::vector<TT> &xValues_){
    if(xValues_.size() != yValues_.size()){
      throw std::logic_error("x and y values list do have the same size.");
    }
    const auto n    = xValues_.size();
    const auto s_x  = std::accumulate(xValues_.begin(), xValues_.end(), 0.0);
    const auto s_y  = std::accumulate(yValues_.begin(), yValues_.end(), 0.0);
    const auto s_xx = std::inner_product(xValues_.begin(), xValues_.end(), xValues_.begin(), 0.0);
    const auto s_xy = std::inner_product(xValues_.begin(), xValues_.end(), yValues_.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
  }
  template <typename T> static double getStdDev(const std::vector<T>& vector_, const std::function<double(const T&)>& evalElementFct_){
    double outVal = 0;
    double mean = getAverage(vector_, evalElementFct_);
    for( auto& element : vector_ ){
      outVal += std::pow(
          static_cast<double>(evalElementFct_(element)) - mean,
          2
      );
    }
    return sqrt( outVal / vector_.size() );
  }


  // Sorting
  template <typename T, typename Lambda> static std::vector<size_t> getSortPermutation(const std::vector<T>& vectorToSort_, const Lambda& firstArgGoesFirstFct_ ){
    std::vector<size_t> p(vectorToSort_.size());
    // 0,1,2,3,4...
    std::iota(p.begin(), p.end(), 0);
    try {
      std::sort(p.begin(), p.end(),
              [&](size_t i, size_t j){
                if( i > vectorToSort_.size() or j > vectorToSort_.size() ) {
                  throw std::runtime_error("std::sort is reaching out of range.");
                }
                return firstArgGoesFirstFct_(vectorToSort_.at(i), vectorToSort_.at(j));
              });
    }
    catch( ... ) {
      // REACHING HERE MIGHT INDICATE THAT SOMETHING IS WRONG WITH firstArgGoesFirstFct_.
      // "not a valid strict weak ordering":
      // Typically: firstArgGoesFirstFct_(obj1, obj2) != !firstArgGoesFirstFct_(obj2, obj1)
      GTLogError << "Something might be wrong with the sort function. Using std::stable_sort instead..." << std::endl;
      try {
        std::stable_sort(p.begin(), p.end(),
              [&](size_t i, size_t j){
                if( i > vectorToSort_.size() or j > vectorToSort_.size() ) {
                  throw std::runtime_error("std::sort is reaching out of range.");
                }
                return firstArgGoesFirstFct_(vectorToSort_.at(i), vectorToSort_.at(j));
              });
      }
      catch( ... ) {
        GTLogError << "Couldn't use std::stable_sort either. Skip sorting." << std::endl;
      }

    }

    return p;
  }
  template <typename T> static std::vector<T> getSortedVector(const std::vector<T>& unsortedVector_, const std::vector<std::size_t>& sortPermutation_ ){
    if(unsortedVector_.empty() or sortPermutation_.size() != unsortedVector_.size()) return {};
    std::vector<T> sortedVec(unsortedVector_.size(), unsortedVector_[0]);
    std::transform(sortPermutation_.begin(), sortPermutation_.end(), sortedVec.begin(),
                   [&](std::size_t i){ return unsortedVector_[i]; });
    return sortedVec;
  }
  template <typename T> static void applyPermutation(std::vector<T>& vectorToPermute_, const std::vector<std::size_t>& sortPermutation_){
    std::vector<bool> done(vectorToPermute_.size(), false);
    for( std::size_t iEntry = 0; iEntry < vectorToPermute_.size(); iEntry++ ){
      if( done[iEntry] ){ continue; }
      done[iEntry] = true;
      std::size_t jPrev = iEntry;
      std::size_t jEntry = sortPermutation_[iEntry];
      while( iEntry != jEntry ){
        std::swap(vectorToPermute_[jPrev], vectorToPermute_[jEntry]);
        done[jEntry] = true;
        jPrev = jEntry;
        jEntry = sortPermutation_[jEntry];
      }
    }
  }
  template <typename T, typename Lambda> static void sortVector(std::vector<T>& vectorToSort_, const Lambda& firstArgGoesFirstFct_){
    std::sort(vectorToSort_.begin(), vectorToSort_.end(), firstArgGoesFirstFct_);
  }
  template <typename T, typename Lambda> static void removeEntryIf(std::vector<T>& vector_, const Lambda& removeIfFct_){
    vector_.erase( std::remove_if(vector_.begin(), vector_.end(), removeIfFct_), vector_.end() );
  }

  // Others
  template<typename T, typename TT> static T& getListEntry(std::list<T>& list_, TT index_){
    typename std::list<T>::iterator it = list_.begin();
    std::advance(it, index_);
    return *it;
  }
  template<typename T, typename TT> static const T& getListEntry(const std::list<T>& list_, TT index_){
    typename std::list<T>::const_iterator it = list_.begin();
    std::advance(it, index_);
    return *it;
  }

  /* ================= CSRVector ================= */

  template <typename T, typename OffsetT>
  CSRVector<T, OffsetT>::CSRVector() = default;

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::reserve(size_t nEntries) {
    offsets_.reserve(nEntries + 1);
  }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::reserve_data(size_t nTotalValues) {
    data_.reserve(nTotalValues);
  }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::clear() {
    offsets_.clear();
    data_.clear();
    openRow_ = false;
  }

  template <typename T, typename OffsetT>
  size_t CSRVector<T, OffsetT>::size() const {
    return offsets_.empty() ? 0 : offsets_.size() - 1;
  }

  template <typename T, typename OffsetT>
  bool CSRVector<T, OffsetT>::empty() const {
    return size() == 0;
  }

  template <typename T, typename OffsetT>
  typename CSRVector<T, OffsetT>::RowBuilder
  CSRVector<T, OffsetT>::emplace_back() {

    if (openRow_)
      throw std::logic_error("CSRVector: only one open row allowed");

    if (offsets_.empty())
      offsets_.push_back(0);

    size_t rowIndex = offsets_.size() - 1;
    offsets_.push_back(static_cast<offset_type>(data_.size()));
    openRow_ = true;

    return RowBuilder(this, rowIndex);
  }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::push_back(const std::vector<value_type>& row) {
    auto rb = emplace_back();
    for (const auto& v : row) rb.emplace_back(v);
    rb.close();
  }

  template <typename T, typename OffsetT>
  typename CSRVector<T, OffsetT>::RowView
  CSRVector<T, OffsetT>::operator[](size_t iEntry) const {

    if (iEntry >= size())
      throw std::out_of_range("CSRVector: row index out of range");

    auto b = static_cast<size_t>(offsets_[iEntry]);
    auto e = static_cast<size_t>(offsets_[iEntry + 1]);
    return RowView(data_.data() + b, e - b);
  }

  template <typename T, typename OffsetT>
  const std::vector<typename CSRVector<T, OffsetT>::offset_type>&
  CSRVector<T, OffsetT>::offsets() const { return offsets_; }

  template <typename T, typename OffsetT>
  const std::vector<typename CSRVector<T, OffsetT>::value_type>&
  CSRVector<T, OffsetT>::data() const { return data_; }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::append(const CSRVector& other) {

    if (openRow_ || other.openRow_)
      throw std::logic_error("CSRVector::append: cannot append with an open row");

    if (other.empty()) return;

    if (this->empty()) {
      offsets_ = other.offsets_;
      data_ = other.data_;
      return;
    }

    const size_t baseData = data_.size();
    const auto base = static_cast<offset_type>(baseData);

    data_.insert(data_.end(), other.data_.begin(), other.data_.end());

    const size_t oldOffsetSize = offsets_.size();
    offsets_.resize(oldOffsetSize + other.size());

    for (size_t i = 0; i < other.size(); ++i) {
      offsets_[oldOffsetSize + i] = base + other.offsets_[i + 1];
    }
  }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::append_move(CSRVector&& other) {

    if (openRow_ || other.openRow_)
      throw std::logic_error("CSRVector::append_move: cannot append with an open row");

    if (other.empty()) return;

    if (this->empty()) {
      offsets_ = std::move(other.offsets_);
      data_ = std::move(other.data_);
      other.openRow_ = false;
      return;
    }

    const size_t baseData = data_.size();
    const auto base = static_cast<offset_type>(baseData);

    // Move elements one by one (works for any T, but still O(n))
    data_.reserve(data_.size() + other.data_.size());
    for (auto& v : other.data_) data_.emplace_back(std::move(v));

    const size_t oldOffsetSize = offsets_.size();
    offsets_.resize(oldOffsetSize + other.size());

    for (size_t i = 0; i < other.size(); ++i) {
      offsets_[oldOffsetSize + i] = base + other.offsets_[i + 1];
    }

    other.clear();
  }

  template <typename T, typename OffsetT>
  CSRVector<T, OffsetT>
  CSRVector<T, OffsetT>::concat(const std::vector<CSRVector>& parts) {

    CSRVector out;

    size_t totalEntries = 0;
    size_t totalData = 0;
    for (const auto& p : parts) {
      if (p.openRow_) throw std::logic_error("CSRVector::concat: part has an open row");
      totalEntries += p.size();
      totalData += p.data_.size();
    }

    out.reserve(totalEntries);
    out.reserve_data(totalData);

    for (const auto& p : parts) out.append(p);

    return out;
  }

  /* ================= RowView ================= */

  template <typename T, typename OffsetT>
  CSRVector<T, OffsetT>::RowView::RowView(const value_type* ptr, size_t n)
    : ptr_(ptr), n_(n) {}

  template <typename T, typename OffsetT>
  const typename CSRVector<T, OffsetT>::value_type&
  CSRVector<T, OffsetT>::RowView::operator[](size_t k) const {
    if (k >= n_)
      throw std::out_of_range("CSRVector: column index out of range");
    return ptr_[k];
  }

  template <typename T, typename OffsetT>
  const typename CSRVector<T, OffsetT>::value_type*
  CSRVector<T, OffsetT>::RowView::data() const { return ptr_; }

  template <typename T, typename OffsetT>
  size_t CSRVector<T, OffsetT>::RowView::size() const { return n_; }

  template <typename T, typename OffsetT>
  bool CSRVector<T, OffsetT>::RowView::empty() const { return n_ == 0; }

  template <typename T, typename OffsetT>
  const typename CSRVector<T, OffsetT>::value_type*
  CSRVector<T, OffsetT>::RowView::begin() const { return ptr_; }

  template <typename T, typename OffsetT>
  const typename CSRVector<T, OffsetT>::value_type*
  CSRVector<T, OffsetT>::RowView::end() const { return ptr_ + n_; }

  /* ================= RowBuilder ================= */

  template <typename T, typename OffsetT>
  CSRVector<T, OffsetT>::RowBuilder::RowBuilder(
    CSRVector* owner, size_t rowIndex)
    : owner_(owner), rowIndex_(rowIndex) {}

  template <typename T, typename OffsetT>
  CSRVector<T, OffsetT>::RowBuilder::RowBuilder(RowBuilder&& other) noexcept {
    owner_ = other.owner_;
    rowIndex_ = other.rowIndex_;
    open_ = other.open_;
    other.owner_ = nullptr;
    other.open_ = false;
  }

  template <typename T, typename OffsetT>
  typename CSRVector<T, OffsetT>::RowBuilder&
  CSRVector<T, OffsetT>::RowBuilder::operator=(RowBuilder&& other) noexcept {
    if (this != &other) {
      close_if_needed();
      owner_ = other.owner_;
      rowIndex_ = other.rowIndex_;
      open_ = other.open_;
      other.owner_ = nullptr;
      other.open_ = false;
    }
    return *this;
  }

  template <typename T, typename OffsetT>
  CSRVector<T, OffsetT>::RowBuilder::~RowBuilder() {
    close_if_needed();
  }

  template <typename T, typename OffsetT>
  template <class... Args>
  typename CSRVector<T, OffsetT>::value_type&
  CSRVector<T, OffsetT>::RowBuilder::emplace_back(Args&&... args) {

    if (!owner_ || !open_)
      throw std::logic_error("CSRVector: cannot emplace_back on closed row");

    owner_->data_.emplace_back(std::forward<Args>(args)...);
    return owner_->data_.back();
  }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::RowBuilder::close() {

    if (!owner_ || !open_)
      throw std::logic_error("CSRVector: row already closed");

    owner_->offsets_[rowIndex_ + 1] =
      static_cast<offset_type>(owner_->data_.size());

    open_ = false;
    owner_->openRow_ = false;
  }

  template <typename T, typename OffsetT>
  void CSRVector<T, OffsetT>::RowBuilder::close_if_needed() {
    if (owner_ && open_) close();
  }


}


#endif // CPP_GENERIC_TOOLBOX_VECTOR_H
