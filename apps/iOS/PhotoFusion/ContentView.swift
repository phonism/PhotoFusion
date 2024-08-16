import SwiftUI
import UIKit
import UniformTypeIdentifiers
import ImageIO
import Photos

struct DocumentPicker: UIViewControllerRepresentable {
    @Binding var selectedData: Data?
    @Environment(\.presentationMode) var presentationMode
    
    class Coordinator: NSObject, UIDocumentPickerDelegate {
        var parent: DocumentPicker
        
        init(parent: DocumentPicker) {
            self.parent = parent
        }
        
        func documentPicker(_ controller: UIDocumentPickerViewController, didPickDocumentsAt urls: [URL]) {
            self.parent.presentationMode.wrappedValue.dismiss()

            if let url = urls.first {
                // 使用 startAccessingSecurityScopedResource 方法
                if url.startAccessingSecurityScopedResource() {
                    defer { url.stopAccessingSecurityScopedResource() }
                    
                    do {
                        let data = try Data(contentsOf: url)
                        self.parent.selectedData = data
                    } catch {
                        print("Error reading file data: \(error)")
                    }
                } else {
                    print("Failed to access security scoped resource")
                }
            } else {
                print("Failed to get file URL")
            }
        }
        
        func documentPickerWasCancelled(_ controller: UIDocumentPickerViewController) {
            self.parent.presentationMode.wrappedValue.dismiss()
        }
    }
    
    func makeCoordinator() -> Coordinator {
        return Coordinator(parent: self)
    }
    
    func makeUIViewController(context: Context) -> UIDocumentPickerViewController {
        let picker = UIDocumentPickerViewController(forOpeningContentTypes: [UTType.image])
        picker.delegate = context.coordinator
        return picker
    }
    
    func updateUIViewController(_ uiViewController: UIDocumentPickerViewController, context: Context) {}
}

struct ZoomableImageView: UIViewRepresentable {
    var image: UIImage
    var dismissAction: () -> Void
    
    func makeUIView(context: Context) -> UIScrollView {
        let scrollView = UIScrollView()
        scrollView.delegate = context.coordinator
        scrollView.maximumZoomScale = 3.0
        scrollView.minimumZoomScale = 1.0
        scrollView.zoomScale = 1.0
        scrollView.showsHorizontalScrollIndicator = false
        scrollView.showsVerticalScrollIndicator = false
        scrollView.contentInsetAdjustmentBehavior = .never
        
        let imageView = UIImageView(image: image)
        imageView.contentMode = .scaleAspectFit
        imageView.isUserInteractionEnabled = true
        scrollView.addSubview(imageView)
        
        context.coordinator.imageView = imageView
        
        // 设置单击手势以关闭视图
        let tapGesture = UITapGestureRecognizer(target: context.coordinator, action: #selector(Coordinator.handleTap(_:)))
        scrollView.addGestureRecognizer(tapGesture)
        
        // 设置imageView的初始化布局约束
        imageView.translatesAutoresizingMaskIntoConstraints = false
        NSLayoutConstraint.activate([
            imageView.widthAnchor.constraint(equalTo: scrollView.widthAnchor),
            imageView.heightAnchor.constraint(equalTo: scrollView.heightAnchor),
            imageView.centerXAnchor.constraint(equalTo: scrollView.centerXAnchor),
            imageView.centerYAnchor.constraint(equalTo: scrollView.centerYAnchor)
        ])
        
        return scrollView
    }
    
    func updateUIView(_ uiView: UIScrollView, context: Context) {
        // 更新视图（如果需要）

    }
    
    func makeCoordinator() -> Coordinator {
        Coordinator(dismissAction: dismissAction)
    }
    
    class Coordinator: NSObject, UIScrollViewDelegate {
        var imageView: UIImageView?
        var dismissAction: () -> Void
        
        init(dismissAction: @escaping () -> Void) {
            self.dismissAction = dismissAction
        }
        
        func viewForZooming(in scrollView: UIScrollView) -> UIView? {
            return imageView
        }
        
        @objc func handleTap(_ gesture: UITapGestureRecognizer) {
            dismissAction()
        }
    }
}

struct ContentView: View {
    @State private var isShowingDocumentPicker = false
    @State private var selectedData: Data? {
        didSet {
            print("selectedData changed")
        }
    }
    @State private var decodedImage: UIImage?
    @State private var outputFileURL: URL? = nil
    @State private var showingSaveSuccessAlert = false
    @State private var isImageViewerPresented = false // 控制 ImageViewer 视图的显示
    
    
    @ViewBuilder
    var body: some View {
        VStack {
            Text("PhotoFusion")
                .font(.largeTitle)
                .padding()
            
            Button(action: {
                self.isShowingDocumentPicker = true
            }) {
                Text("Select RAW File from Files")
                    .padding()
                    .background(Color.blue)
                    .foregroundColor(.white)
                    .cornerRadius(8)
            }
            .sheet(isPresented: $isShowingDocumentPicker) {
                DocumentPicker(selectedData: self.$selectedData)
            }
            
            if selectedData != nil {
                Button(action: {
                    self.decodeRawFile()
                }) {
                    Text("Decode Selected RAW File")
                        .padding()
                        .background(Color.green)
                        .foregroundColor(.white)
                        .cornerRadius(8)
                }
            }
            
            if let image = decodedImage {
                Button(action: {
                    self.isImageViewerPresented = true
                }) {
                    Image(uiImage: image)
                        .resizable()
                        .scaledToFit()
                        .frame(width: 300, height: 300)
                        .border(Color.black, width: 1)
                }
                .fullScreenCover(isPresented: $isImageViewerPresented) {
                    ZoomableImageView(image: image) {
                        // 当用户点击图片时关闭全屏视图
                        self.isImageViewerPresented = false
                    }
                    .edgesIgnoringSafeArea(.all)
                }
                Button(action: {
                    if let fileURL = self.outputFileURL {
                        saveFileToPhotosAlbum(fileURL: fileURL)
                    }
                }) {
                    Text("Save to Photos")
                        .padding()
                        .background(Color.orange)
                        .foregroundColor(.white)
                        .cornerRadius(8)
                }
                .alert(isPresented: $showingSaveSuccessAlert) {
                    Alert(title: Text("Success"), message: Text("Image saved to Photos"), dismissButton: .default(Text("OK")))
                }
            }
        }}
    func saveFileToPhotosAlbum(fileURL: URL) {
        // 请求访问相册权限
        PHPhotoLibrary.requestAuthorization { status in
            guard status == .authorized else {
                print("访问相册的权限未被授权")
                return
            }
            
            PHPhotoLibrary.shared().performChanges({
                // 创建一个照片请求
                let creationRequest = PHAssetCreationRequest.forAsset()
                creationRequest.addResource(with: .photo, fileURL: fileURL, options: nil)
            }) { success, error in
                if success {
                    print("文件已成功保存到相册")
                    DispatchQueue.main.async {
                        self.showingSaveSuccessAlert = true
                    }
                } else {
                    print("保存到相册失败: \(String(describing: error))")
                }
            }
        }
    }
    
    func decodeRawFile() {
        guard let data = selectedData else { return }
        do {
            // 将数据写入到临时文件
            let tempURL = FileManager.default.temporaryDirectory.appendingPathComponent("temp.nef")
            try data.write(to: tempURL)
            
            // 使用临时文件解码
            // try RDLibRawDecoder.decodeRawImage(fromFileURL: tempURL)
            
            // 获取文档目录路径
            let tempDirectoryURL = tempURL.deletingLastPathComponent()
            let outputFileURL = tempDirectoryURL.appendingPathComponent("output.tiff")
            let raw2TiffWrapper = Raw2TiffWrapper()

            raw2TiffWrapper.convertRaw(toTiff: tempURL.path, tiffFilePath: outputFileURL.path)
            self.outputFileURL = outputFileURL // 更新文件路径
            // 读取保存的 TIFF 文件并显示
            let tiffData = try Data(contentsOf: outputFileURL)
            
            // 使用 UIImage 直接加载 TIFF 数据
            if let uiImage = UIImage(data: tiffData) {
                self.decodedImage = uiImage
            } else {
                print("Failed to create UIImage from TIFF data")
            }
        } catch {
            print("Failed to decode image: \(error)")
        }
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: ContentView {
        ContentView()
    }
}
