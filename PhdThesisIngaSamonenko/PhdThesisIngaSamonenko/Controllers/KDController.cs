using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class KDController : Controller
    {
        public ActionResult Index()
        {
            return File("KD.R", "text/plain");
        }

    }
}
